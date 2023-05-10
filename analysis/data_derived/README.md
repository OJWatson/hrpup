## File/Data Description

### HRP2 Risk Map Objects

1. `R6_map.rds` - R6 map object for plotting risk maps. See `?hrpup:::R6_hrp2_map` for details

This object can be used to produce risk maps for specific scenarios for specific global regions and for either the innate or composite (Africa only) risk score:

```
hrp2_map$plot(
  region = "africa", 
  Micro.2.10 = "central", 
  ft = "central", 
  microscopy.use = "central", 
  rdt.nonadherence = "central", 
  fitness = "central", 
  rdt.det = "central",
  risk = "innate"
)
```

### Parameter Files for Simulations Specific to DIDE Cluster

1. `param_start.rds` - Parameters for burnin runs
1. `param_grid.rds` - Parameters for grid continuations
1. `final_param_grid.rds` - Parameters for grid continuations with IDs of burnin runs to pick up
