## Impact Analysis Data OutputsFile/Data Description

This directory, `data_out`, within `impact_analysis` directory contains outputs 
from modelling the public health impact of the spread of **pfhrp2/3 deletions** 
across sub-Saharan Africa. These results estimate changes in malaria burden 
associated with continued use of HRP2-based rapid diagnostic tests (RDTs) 
compared to scenarios where countries transition to combined LDH/HRP2-based RDTs 
upon reaching the WHO-defined 5% false-negative threshold.

### File Structure

All `.csv` files in this directory share a common structure with columns detailed below.

### Identifiers
- **`id_1`**: Unique geographic identifier for each first-level administrative region.
- **`t`**: Simulated time in years from from t = 0 (2023) to t = 20 years. 
- **`delay`**: Number of years delay between reaching the 5% threshold and switching to combined LDH/HRP2-based RDTs.
- **`scenario`**: Different scenarios of the speed of HRP2 selection due to uncertainty in parameter estimates of risk factors of selection (e.g. malaria prevalence, treatment seeking rates etc).
- **`type`**: Type of policy/intervention scenario: `"5% Threshold Strategy"` or `"No RDT Switching"` or `"Switch Now"`. The `"5% Threshold Strategy"` models admin regions swapping to alternative RDTs when the 5% false-negative RDTs threshold is crossed with set delays. The `"No RDT Switching"` models a counterfactual in which countries continue to use current RDTs. The `"Switch Now"` models a scenario in which everywhere switched at t = 0.  

### Malaria Burden Metrics

Each malaria burden metric includes lower (lci), median (med), and upper (uci) estimates, representing uncertainty intervals (2.5%, 50%, 97.5%) across stochastic simulations.

| Column Name               | Description                                             |
|---------------------------|---------------------------------------------------------|
| **`freq_lci/med/uci`**            | False-negative RDT frequency estimates                  |
| **`micro_2_10_lci/med/uci`**      | Malaria prevalence (microscopy detectable, ages 2–10)   |
| **`pcr_lci/med/uci`**             | Malaria prevalence (PCR detectable, all ages)           |
| **`clinical_05_lci/med/uci`**     | Clinical malaria incidence (children under 5)           |
| **`clinical_lci/med/uci`**        | Clinical malaria incidence (all ages)                   |
| **`severe_05_lci/med/uci`**       | Severe malaria cases (children under 5)                 |
| **`severe_lci/med/uci`**          | Severe malaria cases (all ages)                         |
| **`mortality_05_lci/med/uci`**    | Malaria mortality (children under 5)                    |
| **`mortality_100_lci/med/uci`**   | Malaria mortality (all ages)                            |

## Extra Information

- In all the outputs above, the `scenario` identifier covers three scenarios. The **central** times provide the estimated times based on the central parameter estimate for each of the parameters that we explored and that are known to impact the speed of selection of ArtR (malaria prevalence, treatment seeking related parameters, diagnostic usage and adherence, cross-reactivity of HRP3, fitness costs). The **best** times assume the upper or lower estimate value for each parameter (depending on the direction of its effect on selection) such that the selection of ArtR will increase at its slowest. Conversely, the **worst** times assume the upper or lower value for each parameter (depending on the direction of its effect on selection) such that the selection of ArtR will increase at its fastest.
- All numeric outputs have been rounded to 6 significant figures to reduce memory of the files exported.
- Each output file has the `"No RDT Switching"` scenario (which has `delay` labelled as 100 and identified by `type`), the `"Switch Now"` scenario (`delay` labelled as 0 and identified by `type`) and the delay scenario in the `"5% Threshold Strategy"`.
- The `<>_wmrscaled_<>` files have been calibrated and scaled to align with WMR estimated malaria case incidence and mortality
