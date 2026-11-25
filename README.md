# Train-bridge analysis

This script provides a simple way to analyse the bending moments along a bridge's surface from a train's load, modeled as point loads at the wheels.

<u>Current features</u>:

- Bending moment diagram (BMD) for given train dimensions, weight, bridge length, and train position.
- 3D graph comparing BMD to train position, to analyze maximum possible bending moment.
- Cross sectional analysis across non uniform bridge including I, centroid
- Analysis of maximum bending stress for all train positions
- Analysis of maximum shear stress for all train positions

_note: **get_constants.m** contains sample data, though you can easily replace it to match your situation as long as everything is consistent._
