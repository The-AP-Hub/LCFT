# LCFT
The L'Varian Coherence Field Theory

This repository contains basic stubs for the coherence-field implementation. The `src/` directory provides
simple C++ classes:

- `Point` – a lightweight 2D point structure shared across the codebase.
- `CoherenceField` – represents a coherence field and exposes a method to sample the field at a point.
- `AdaptiveMesh` – maintains a set of mesh nodes and allows refinement around a given point.
- `SystemMatrices` – stores matrices used by the simulation in a flattened form.

These files serve as a starting point for further development.
