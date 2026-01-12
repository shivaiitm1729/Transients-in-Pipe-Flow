Overview

This repository contains a numerical solver for unsteady (transient) flow in pressurized pipe systems, commonly known as water hammer. The model simulates pressure and velocity transients caused by gradual valve closure in a reservoir–pipe–valve system.

The problem setup and parameters are based on the classical example presented in Applied Hydraulic Transients by M.H.Chaudhry (page 104), including:

Two pipes with different geometric and hydraulic properties. A downstream valve with a prescribed closure law and a computation of transient head and discharge variations

Problem Description:

The system consists of an upstream reservoir with constant head Pipe No. 1 connected to Pipe No. 2. A downstream valve undergoing controlled closure with an initial steady-state flow prior to valve motion

Pipe Properties:

Pipe No. 1

Length: 550 m

Diameter: 0.75 m

Wave speed: 1100 m/s

Darcy friction factor: 0.010

Pipe No. 2

Length: 450 m

Diameter: 0.60 m

Wave speed: 900 m/s

Darcy friction factor: 0.012

Description of valve closure

Valve closure follows a prescribed closure curve over approximately 6 seconds

Closure law is smooth and time-dependent, avoiding instantaneous shutoff

Valve motion directly controls discharge and induces pressure waves
