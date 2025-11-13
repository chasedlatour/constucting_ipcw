## Constructing inverse probability of censoring weights (IPCW) with R

We supply three scripts to demonstrate an example of IPCW in R:

-   `ipcw_example.R` is the primary outward facing script that calls the functions in the others.

-   `helper_functions.R` contains functions to help with the data preparation for the example:

    -   `cutoff_followup`: a function for reading in the data and re-deriving censoring and event indicators based on a user-specified cut off time.

    -   `qrspline`: a function to derived quadratic restricted splines to allow for flexibly modeling time and continuous covariates.

-   `ipcw_functions.R` contains function to assist in stepping through the IPCW process:

    -   `convert_to_long`: converts a conventional time-to-event dataset (i.e., one row per observational unit with a variable *t* and event indicator *delta*) to a long data frame with a row per time interval per observational unit. Users can change the interval length (default 1 day) and specify the censoring and event indicators.

    -   `ipcw`: takes in a long data frame with discrete time intervals and time-varying event and censoring indicators and a user-specified model form and estimates the unstabilized inverse probability of censoring weights.

    -   `compare_cumul_inc`: estimates the IPC weighted and naive survival functions and returns a list with an element with the cumulative incidence functions for each and a simple plot comparing the curves. Requires both the original conventional time-to-event data frame and the data frame in long format.
