# FEMA_P58_damage_correlation

This simple Matlab code demonstrates a proposed algorithm for producing simulations of partially dependent component damage in the FEMA P-58 performance assessment methodology.

The algorithm demonstrated here is described in the following document

> Jack W. Baker, Ed Almeter, Dustin Cook, Abbie Liel, and Curt Haselton (2023) "A model for partially dependent component damage fragilities in seismic risk analysis," (in review).

The results from this calculation are the "Simple Example" presented
in the above manuscript.

## Organization

'MAIN_correlations.m' is a script where the input parameters for the simulation are specified. The script correctly formats the parameters and then calls the analysis function.

'fn_simulate_damage.m' is the analysis function that performs the component damage simulations. Component damage realizations are produced, summary statistics are computed, and an annotated figure is produced to summarize the results.
