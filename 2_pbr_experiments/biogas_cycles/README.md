# Biogas cycling experiment

Data and analysis for three countercurrent synthetic-biogas cycles on 31 May 2023, using a ceria packed bed at 930 °C.

| File | Description |
| --- | --- |
| `rx_test_027.lvm` | Original instrument data (tab-separated, decimal commas). |
| `plot_data.py` | Extracts measurements, corrects the H2 signal, and saves the processed CSV and plots. |
| `31052023_biogas_930C_dataset.csv` | Processed time series of temperature, pressure, gas flows, concentrations, and flow set points. |
| `analysis_biogas.py` | Integrates cycle flows and calculates conversions, yields, elemental balances, and model feed inputs. |
| `analysis_biogas_output.txt` | Saved analysis output, including per-cycle values and three-cycle totals. |
| `analysis_biogas_summary.csv` | Summary metrics; gas inlet amounts are three-cycle totals. |
| `reactor_bed` | LaTeX table of bed dimensions, ceria loading, operating conditions, and void fractions. |
| `plots/31052023_Biogas_930C_dataset.png` | Overview plot of gas concentrations, inlet flows, and bed temperature. |
| `plots/31052023_Biogas_930C_dataset.pdf` | PDF version of the overview plot. |
