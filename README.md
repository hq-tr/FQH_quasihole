# FQH_quasihole
(Beta version. Last updated Feb 2022)

Main function: `FQH_quasiholes.FQH_quasihole_poly(Ne, pos, ground_state)`

Add fluxes at positions pos=[w_1, w_2,...] into a given many_body state.

Function arguments:

- `Ne`: Number of electrons
- `pos`:list of complex numbers that indicate positions of the added flux 
- `ground_state`: an `FQH_states.fqh_state` variable

This functions requires `FQH_states.py` and `misc.py` from _qhe-library.


