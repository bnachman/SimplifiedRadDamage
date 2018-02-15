
## Simplified Rad Sim

This package employs a simplified model of silicon radiation damage to estimate the charge collection efficiency inside a planar sensor.  The code is all contained in the file is RadHelpers.py and two example programs called Example_IntegratedFluence and Exapmle_PerFluence are provided.  RadHelpers does not depend on ROOT but the examples depend on ROOT for plotting.  

## Physics

This simplified model assumes a uniform electric field inside the sensor bulk and models the Ramo potential as a double-exponential, which is known to be an excellent approximation based on analytic and numerical calculations.  Charge trapping is modeled with a single rate for electrons and for holes.  Since charge trapping and charge induction are described by exponentials, it is simple to analytically describe the entire model for the collected electron and hole charge as a function of depth inside the sensor.  The induced charge on neighboring pixels as well as any transverse motion of charge carriers is ignored in the model.