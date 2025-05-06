WISDOM
    Weighted
    Ion
    Simulations of
    Diagnostic
    Orbit
    Space
    Maps

Run instructions:
    NOT CURRENTLY FUNCTIONAL IN THIS STATE - USER BEWARE

Contents:

    ------------------------------------------------
    Julia code
    ------------------------------------------------

    atomic.jl     : Atomic cross section data
    beams.jl      : Neutral beam formalism
    collisions.jl : Collision operators
    IO.jl         : input/output functions.
    ions.jl       : Ion functions, includes MEP & Boris Push
    main.jl       : Main file for running code. Called by 
    plots.jl      : Plotting functions
    setup.jl      : Set up script (run on each worker).
    sigmav.jl     : Functionality for pre-computing distribution-averaged sigma-v tables. **NOT IMPLEMENTED**
    utils.jl      : Utility functions - equilibrium models, math, and array operations.
    weights.jl    : Functions for building single particle weights for each detector.

    ------------------------------------------------
    Other files
    ------------------------------------------------

    run_local.sh  : Bash script, creates and evaluates a job locally.
    *.toml        : Julia project files
 
Contact: Violet Player, violetmplayer@gmail.com
