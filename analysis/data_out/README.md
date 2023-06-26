## File/Data Description

The `data_out` directory contains modelled outputs for external partners (as opposed to results directly used for
risk mapping of *pfhrp2/3* deletions) at the first administrative unit. 

### Modelled risk outputs

The first set of outputs provide the estimated time that false-negative HRP2-based RDTs due to *pfhrp2/3* deletions will 
increase from 1% to 5%. Where times are predicted to take longer than 40 years, these have been censored (">40"):

1. `central_times.csv`
1. `optimistic_times.csv`
1. `pessimistic_times.csv`

The second set of outputs provide the estimated time from now that false-negative HRP2-based RDTs due to *pfhrp2/3* deletions will 
reach 5% based on the current best estimates of the prevalence of the deletions.

These time are only provided for African regions because this is a) where the main malaria burden and testing via HRP2-based RDTs occurs, but also where the greatest resolution on *pfhrp2/3* prevalence is available due to the recent increase in conducted surveys. 

Similarly to before, where times are predicted to take longer than 40 years, these have been censored (">40"):

1. `central_times_prospective.csv`
1. `optimistic_times_prospective.csv`
1. `pessimistic_times_prospective.csv`

In all the outputs above, the **central** times provide the estimated times based on the central parameter estimate for each of the parameters that we explored and that are known to impact the speed of selection for *pfhrp2/3* deletions (malaria prevalence, treatment related parameters, fitness costs and HRP3 cross-reactivity factors). The **optimistic** times assume the 2.5% or 97.5% percentile value for each parameter (depending on the direction of its effect on selection) such that the selection of *pfhrp2/3* deletions will increase at its slowest. Conversely, the **pessimistic** times  assume the 2.5% or 97.5% percentile value for each parameter (depending on the direction of its effect on selection) such that the selection of *pfhrp2/3* deletions will increase at its fastest.

Lastly, the uncertainty provided in each of the files, reflects the uncertainty that arises from our use of a stochastic, individual based model for the selection of *pfhrp2/3* deletions. 

### Overall risk classifications

1. `full_results.csv`: Provides the discretised risk category for each administrative unit across each our parameter combinations. 

