Phase I — Numerical / Plasma Core
────────────────────────────────────
M1  Mesh / Field / Boundary              ✅
M2  Linear algebra                       ✅
M3  Poisson                              ✅
M4  Diffusion / advection                ✅
M5  Scharfetter-Gummel                   ✅
M6  Species continuity                   ✅
M7  Electrostatic coupling               ✅
M8  Multi-species + ReactionNetwork      ✅
M9  Simulation / TimeLoop                ✅


Phase II — Robust Simulation
────────────────────────────────────
M10 Adaptive timestep                   ✅
M11 Diagnostics / Observability         ✅
    M11.1 canonical unit convention - use mp-units ✅
    M11.2 2D physical-volume semantics ✅
    M11.3 scalar field statistics ✅
    M11.4 species diagnostics ✅
    M11.5 charge diagnostics ✅
    M11.6 field diagnostics ✅
    M11.7 timestep diagnostics ✅
    M11.8 structured logging ✅
    M11.9 failure diagnostics ✅

M12 Output / checkpoint                 ✅


Phase III — Better Plasma Physics
────────────────────────────────────
M13 Electron energy equation            ✅
    M13.1 energy quantities / metadata  ✅
    M13.2 mean-energy conversion        ✅
    M13.3 SG energy transport           ✅
    M13.4 explicit CFL / positivity     ✅
    M13.5 numerical tests / document    ✅
    M13.6 Simulation integration        ✅
M14 E/N / Te dependent chemistry        ✅
    M14.1 reduced-field / Te conversion ✅
    M14.2 tabulated rate coefficients   ✅
    M14.3 cell-local reaction-rate assembly ✅
    M14.4 numerical tests / document    ✅
    M14.5 Simulation integration        ✅
M15 Plasma wall / electrode BC             ✅
    M15.1 wall particle/energy flux contract ✅
    M15.2 secondary-emission flux assembly ✅
    M15.3 continuity-equation wall discretization ✅
    M15.4 electron-energy wall discretization ✅
    M15.5 numerical tests / document      ✅
    M15.6 Simulation integration          ✅


Phase IV — Validation
────────────────────────────────────
M16 Plasma benchmark suite


Phase V — Performance
────────────────────────────────────
M17 Backend abstraction refinement
M18 GPU implementation
M19 GPU/CPU numerical equivalence
M20 Performance optimization
