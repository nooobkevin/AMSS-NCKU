# AMSS-NCKU Repository Analysis

## Repository Overview

**AMSS-NCKU** is a numerical relativity program developed in China for simulating **binary black hole (BBH) mergers** and computing gravitational wave signatures. It uses the **finite difference method** with **adaptive mesh refinement (AMR)** to numerically solve Einstein's equations.

---

## Phase 1: What It Does (Inputs and Interface)

### Program Purpose
- Solves Einstein Field Equations for binary/multiple black hole systems
- Calculates time evolution of spacetime geometry
- Computes gravitational wave emission during merger events
- Tracks black hole trajectories and horizon properties

### Key Input Parameters

#### Physical Parameters (in `AMSS_NCKU_Input.py`)
- **Black Hole Properties**: 
  - Mass ratio (auto-rescaled to total M=1)
  - Dimensionless spin vectors (aₓ, aᵧ, aᵧ) per BH
  - Initial positions (Brugmann convention: y-axis aligned)
  - Initial momenta (3D vectors)
  - Optional: charge for BSSN-EM variant

- **System Configuration**:
  - Number of punctures (typically 2 for BBH)
  - Initial orbital distance and eccentricity
  - Symmetry: equatorial/octant/none

#### Numerical Parameters
- **Grid Structure**:
  - 9 AMR levels (5 static, 4 moving)
  - 96-120 grid points per direction
  - Largest box: 320×320×320 units
  - Resolution ratio: 2:1 between levels
  - Grid types: "Patch" or "Shell-Patch"

- **Time Evolution**:
  - Evolution time: 0-200 M
  - Courant factor: 0.4
  - Finite difference order: 2nd/4th/6th/8th
  - Time integrator: 4th-order Runge-Kutta
  - Dissipation: Kreiss-Oliger (ε=0.15)

- **Equation Choice**:
  - BSSN (vacuum general relativity)
  - Z4C (constraint-damping variant)
  - BSSN-EM (electrovacuum)
  - BSSN-EScalar (scalar-tensor / f(R) theories)

- **MPI/GPU**:
  - Configurable MPI processes (default: 8)
  - Optional GPU acceleration (currently buggy)

### Workflow (`AMSS_NCKU_Program.py`)

1. **Setup**: Parse input → Generate directories → Validate parameters
2. **Grid Generation**: Create AMR hierarchy → Plot grid structure
3. **Macro Generation**: Auto-generate `macrodef.h` and `macrodef.fh` based on finite-difference order
4. **Compilation**: Copy source → Compile `ABE` (main code) and `TwoPunctureABE` (initial data)
5. **Initial Data**: Run TwoPuncture to solve initial value problem
6. **Time Evolution**: Run ABE to evolve Einstein equations
7. **Output**: Binary dumps, GW waveforms, constraint violations
8. **Post-processing**: Plot orbits, waveforms, ADM mass, constraints

### Output Products
- Black hole trajectories (`bssn_BH.dat`)
- Gravitational wave Ψ₄ at 11 detector locations (50-150 M)
- ADM mass and angular momentum (`bssn_ADMQs.dat`)
- Constraint violations (`bssn_constraint.dat`)
- 2D slices and 3D binary data dumps

---

## Phase 2: Entry Point and Main Loop (Skeleton)

### Entry Point: `ABE.C`

**main() Function** (lines 63-477):
1. **Initialize MPI**: Set up parallel communication
2. **Parse Parameters**: Read `input.par` (or command-line specified file)
3. **Create Evolution Object**: Instantiate `bssn_class` (or variant based on `ABEtype`)
   ```cpp
   #if (ABEtype == 0)
   ADM = new bssn_class(Courant, StartTime, TotalTime, ...);
   #elif (ABEtype == 1)  // BSSN-EScalar
   #elif (ABEtype == 2)  // Z4C
   #elif (ABEtype == 3)  // BSSN-EM
   ```
4. **Setup Initial Data**:
   - `Setup_Initial_Data_Cao()` - Analytical (Bowen-York)
   - `Setup_KerrSchild()` - Kerr-Schild analytical
   - `Read_Ansorg()` - Numerical (TwoPuncture output)
5. **Main Evolution**: `ADM->Evolve(Steps)`
6. **Cleanup**: Finalize MPI and deallocate

### Core Data Structure: `bssn_class.h`

**Key Member Variables**:
- **Evolution State**:
  - `PhysTime` - current physical time
  - `dT` - base timestep (Courant-limited)
  - `Steps` - iteration counter

- **Grid Hierarchy**:
  - `cgh *GH` - Computational Grid Hierarchy
  - `ShellPatch *SH` - Shell coordinate patches

- **BSSN Fields** (3 time levels: current, previous, RHS):
  - `phi` (conformal factor χ = e^(-4φ))
  - `trK` (trace of extrinsic curvature K)
  - `gxx, gxy, gxz, gyy, gyz, gzz` (conformal 3-metric γ̃ᵢⱼ)
  - `Axx, Axy, Axz, Ayy, Ayz, Azz` (traceless extrinsic curvature Ãᵢⱼ)
  - `Gmx, Gmy, Gmz` (auxiliary connection functions: Gmx=Γ̃ˣ, Gmy=Γ̃ʸ, Gmz=Γ̃ᶻ)
  - `Lap` (lapse function α)
  - `Sfx, Sfy, Sfz` (shift vector: Sfx=βˣ, Sfy=βʸ, Sfz=βᶻ)
  - `dtSfx, dtSfy, dtSfz` (shift drivers ∂ₜβⁱ)

- **Diagnostic Fields**:
  - Hamiltonian constraint
  - Momentum constraints
  - Ricci tensor components

### Main Evolution Loop: `bssn_class.C::Evolve()`

**Structure** (lines 2001-2268):
```
for (step = 0; step < MaxSteps; step++):
    ├─ RecursiveStep(level=0)              // Multi-level time stepping
    │  ├─ For coarse level:
    │  │  └─ Step(lev, YN)                 // Single level RK4 step
    │  │     ├─ f_rungekutta4_rout()       // 4-stage RK4 integration
    │  │     │  └─ compute_rhs_bssn()      // Compute ∂ₜq = F(q)
    │  │     ├─ f_sommerfeld_routbam()     // Boundary conditions
    │  │     └─ ghost_zone exchange         // MPI communication
    │  ├─ For fine levels (recursive):
    │  │  ├─ Multiple substeps (2×)
    │  │  ├─ Prolongation (coarse→fine)
    │  │  └─ Restriction (fine→coarse)
    │  └─ RestrictProlong()                // Sync AMR levels
    │
    ├─ PhysTime += dT                      // Advance time
    │
    ├─ AnalysisStuff()                     // Periodic analysis
    │  ├─ Compute constraints
    │  ├─ Extract Ψ₄ (gravitational waves)
    │  ├─ Track black hole positions
    │  └─ Surface integrals (ADM mass)
    │
    ├─ Data Output (if time thresholds met):
    │  ├─ Dump binary data
    │  ├─ Write 2D slices
    │  └─ Checkpoint for restart
    │
    └─ Check termination conditions
```

**RK4 Integration Details** (lines 3606-3610):
```fortran
call f_rungekutta4_rout(
    shape,                   ! grid dimensions
    dT_lev,                 ! timestep for this level
    q_current,              ! current state vector
    q_intermediate,         ! RK4 scratch space
    q_rhs,                  ! right-hand side F(q)
    iter_count)             ! RK substage (1-4)
```

**Runge-Kutta Substages**:
1. **k₁**: Evaluate RHS at t_n, q_n
2. **k₂**: Evaluate RHS at t_n + dt/2, q_n + dt/2 × k₁
3. **k₃**: Evaluate RHS at t_n + dt/2, q_n + dt/2 × k₂
4. **k₄**: Evaluate RHS at t_n + dt, q_n + dt × k₃
5. **Update**: q_{n+1} = q_n + dt/6 × (k₁ + 2k₂ + 2k₃ + k₄)

---

## Phase 3: Physics Modules (Implementation)

### Initial Data: `TwoPunctures.C`

**Purpose**: Solve the Einstein constraint equations to generate initial data for BBH systems using the **puncture method** (Ansorg et al. 2004).

**Method**:
1. **Conformal Decomposition**: Write metric as ψ⁴ times a flat background
2. **Bowen-York Extrinsic Curvature**: Encode linear/angular momentum
3. **Spectral Elliptic Solve**: Find conformal factor ψ satisfying Hamiltonian constraint
4. **Iterate**: Adjust bare masses to match target ADM masses

**Key Variables Set**:
- Initial conformal factor φ
- Initial metric γ̃ᵢⱼ
- Initial extrinsic curvature Kᵢⱼ
- Puncture locations and parameters

### Time Evolution: `bssn_rhs.f90`

**Purpose**: Compute ∂ₜq = F(q) for the 17 BSSN evolution variables.

**BSSN Formulation** (Baumgarte-Shapiro-Shibata-Nakamura):
- Evolves conformal 3-metric γ̃ᵢⱼ (det=1) instead of physical metric
- Introduces auxiliary variable Γ̃ⁱ to make system strongly hyperbolic
- Uses conformal rescaling: χ = e^(-4φ), γ̃ᵢⱼ = χ γᵢⱼ

**Evolution Equations** (abbreviated):
```
∂ₜφ = -1/6 α K + β^i ∂ᵢφ
∂ₜK = -γ^ij D_i D_j α + α(Ãᵢⱼ Ãⁱʲ + 1/3 K²) + 4π α(ρ + S)
∂ₜγ̃ᵢⱼ = -2α Ãᵢⱼ + β^k ∂_k γ̃ᵢⱼ + ...
∂ₜÃᵢⱼ = χ[−DᵢDⱼα + α(Rᵢⱼ − 8πSᵢⱼ)]^TF + β^k ∂_k Ãᵢⱼ + ...
∂ₜΓ̃ⁱ = −2Ãⁱʲ∂ⱼα + 2α(Γ̃ⁱⱼₖ Ãʲᵏ − 2/3 γ̃ⁱʲ∂ⱼK − 8π γ̃ⁱʲSⱼ) + ...
```

**Gauge Conditions** (8 variants supported):
- **Lapse**: 1+log slicing, harmonic, maximal
- **Shift**: Gamma-driver (moving puncture), harmonic

**Key Computations in `compute_rhs_bssn()`** (~1700 lines):
1. **Inverse metric**: gup^ij from g_ij (lines 180-201)
2. **First derivatives**: ∂ᵢ(all fields) via finite differences
3. **Christoffel symbols**: Γⁱⱼₖ from metric derivatives (lines 235-257)
4. **Ricci tensor**: R̃ᵢⱼ from second derivatives (lines 377-443)
5. **Extrinsic curvature**: Raise indices on Ãᵢⱼ (lines 260-282)
6. **RHS assembly**: Combine all terms per evolution equation

**Problem Identified**: **Massive code repetition** in tensor operations:
- Lines 205-231: Manual tensor contractions for Γ̃ⁱ residual (27 terms each, 3 components)
- Lines 260-282: Raising indices A^ij = g^ik g^jl A_kl (6 components, ~8 terms each)
- Lines 326-353: Contractions Γ^i_jk with various tensors
- Lines 356-375: First Christoffel symbols from metric
- Lines 379-400: Ricci tensor computation

**Potential for Optimization**: All these operations are **tensor contractions** that could be expressed concisely using `numpy.einsum()` or similar tools, reducing ~500 lines of error-prone manual indexing to ~20 lines.

### Time Integrator: `rungekutta4_rout.f90`

**Purpose**: Generic 4th-order Runge-Kutta integrator.

**Implementation**:
```fortran
subroutine f_rungekutta4_rout(ex, dT, f0, f1, rhs, YN)
  ! YN=1: f1 = f0 + 0.5*dT*rhs     (k1 → k2)
  ! YN=2: f1 = f0 + 0.5*dT*rhs     (k1 → k3)
  ! YN=3: f1 = f0 + 1.0*dT*rhs     (k1 → k4)
  ! YN=4: f1 = f0 + dT/6*(k1 + 2k2 + 2k3 + k4)
end subroutine
```

**Features**:
- Vectorized over full 3D grid arrays
- Supports scalar, complex, and array-valued fields
- Called once per RK substage for each variable

---

## Phase 4: Infrastructure (Grid Management)

### Computational Grid Hierarchy: `cgh.C/h`

**Purpose**: Manage multi-level adaptive mesh refinement.

**Key Components**:
- **Grid Levels**: Array of `grids[0...maxl-1]` storing number of patches per level
- **Bounding Boxes**: `bbox[level][patch]` - spatial extent of each patch
- **Patch Lists**: `PatL[level]` - linked list of all patches on that level
- **Resolution**: Each finer level has 2× resolution of parent

**Key Methods**:
- `Regrid()`: Dynamically move refined grids to follow black holes
- `recompose_cgh()`: Rebuild grid hierarchy after regridding
- `SetupGhostZones()`: Configure inter-patch communication
- `Prolongation/Restriction`: Transfer data between levels

**Parallelization**: 
- Patches distributed across MPI ranks
- Ghost zone exchanges via MPI communication
- Multiple strategies: domain decomposition, level-based

### Patch Structure: `patch.C/h, MPatch.C`

**Patch Class**:
```cpp
class patch {
    int lev;                    // AMR level
    int shape[3];               // (nx, ny, nz)
    double bbox[6];             // (xmin, xmax, ymin, ymax, zmin, zmax)
    double ***data;             // 3D arrays for variables
    Block *blb, *ble;           // Linked list of blocks
    ...
};
```

**Block Structure**: 
- Subdivides patch into smaller domains for load balancing
- Each block assigned to an MPI rank
- Contains local data arrays and buffer zones

**Key Operations**:
- `Interp_Points()`: Interpolate from coarse to fine grid
- `Restrict()`: Average from fine to coarse grid
- `ApplyStencil()`: Apply finite difference operators
- `FillGhostZones()`: Exchange boundary data

### Shell Patches: `ShellPatch.C/h`

**Purpose**: Handle wave extraction zones in spherical coordinates.

**Motivation**: 
- Far-field gravitational waves best computed on spheres
- Inner region uses Cartesian patches
- Shell patches transition between coordinate systems

**Structure**:
- **Spherical Shells**: r ∈ [r_inner, r_outer], θ ∈ [0,π], φ ∈ [0,2π]
- **Sub-patches**: 1-6 patches (±x, ±y, ±z) depending on symmetry
- **Overlap Regions**: Exchange data with Cartesian patches

**Key Methods**:
- `Coordinates_Cart2Shell()`: Cartesian → (r, θ, φ) transformation
- `Jacobian_Shell2Cart()`: Coordinate transformation Jacobians
- `Interp_Shell_Cartesian()`: Interpolate between coordinate systems
- `Extract_Psi4()`: Compute Newman-Penrose scalar on sphere

**Gravitational Wave Extraction**:
- Compute Weyl curvature Ψ₄ on each shell
- Decompose into spin-weighted spherical harmonics
- Output ₂Y_lm coefficients for post-processing

---

## Key Findings and Problems

### 1. **Code Quality Issues** (Legacy Codebase)
- Written in 2007, uses outdated C++ (pre-C++11, lacks auto, lambda, smart pointers, move semantics)
- Would benefit from modern C++17/C++20 features (structured bindings, ranges, concepts)
- Mixed C++/Fortran90 with manual interfaces (could use ISO_C_BINDING)
- Minimal comments, cryptic variable names
- No modern build system (handwritten Makefiles, would benefit from CMake)
- High technical debt accumulated over ~18 years

### 2. **Massive Code Duplication**
- `bssn_rhs.f90` (1700 lines), `Z4c_rhs.f90` (1800 lines), `bssnEScalar_rhs.f90` (1700 lines)
- Nearly identical structure, differing only in equation terms
- Tensor contractions written out explicitly (~500 lines per file)
- Same pattern repeated in every equation variant

### 3. **Manual Tensor Algebra** 
Example from `bssn_rhs.f90` (lines 260-282):
```fortran
! Raise indices: A^ij = gup^ik gup^jl A_kl
Rxx = gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)
Ryy = gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)
! ... 4 more components
```

**Better approach** (Python/NumPy):
```python
# Define metric as 3×3 array
g_up = np.array([[gupxx, gupxy, gupxz],
                 [gupxy, gupyy, gupyz],
                 [gupxz, gupyz, gupzz]])
A_down = np.array([[Axx, Axy, Axz],
                   [Axy, Ayy, Ayz],
                   [Axz, Ayz, Azz]])

# One-liner: A^ij = gup^ik gup^jl A_kl
A_up = np.einsum('ik,jl,kl->ij', g_up, g_up, A_down)
```

### 4. **Parallel Infrastructure Complexity**
- Custom AMR implementation (not using established libraries like Carpet, GRChombo)
- Manual MPI ghost zone exchanges
- Complex patch decomposition logic
- Difficult to debug and extend

### 5. **GPU Support Issues**
- GPU code (`bssn_gpu.cu`) is "buggy" per comments
- Kernel launches hardcoded for specific grid sizes
- No performance comparisons with CPU code

---

## Refactoring Recommendations

### Short-term (Maintainability)

1. **Consolidate RHS Functions**
   - Create generic `tensor_contraction()` utilities
   - Use Fortran intrinsic `MATMUL` or external BLAS
   - Reduce ~1500 lines to ~300 per RHS file

2. **Add Documentation**
   - Inline comments explaining equations
   - Map variable names to physics notation
   - Add references to papers

3. **Modernize Build System**
   - Use CMake instead of handwritten Makefiles
   - Automate macro generation
   - Better compiler flag management

### Mid-term (Performance)

4. **Adopt Standard AMR Framework**
   - Consider porting to Carpet (Cactus)
   - Or use GRChombo (Chombo AMR + BSSN)
   - Leverage community-tested infrastructure

5. **Fix GPU Implementation**
   - Profile GPU vs CPU performance
   - Use portable GPU frameworks (Kokkos, RAJA)
   - Auto-generate kernels from equations

6. **Optimize Tensor Operations**
   - Use optimized linear algebra (BLAS/LAPACK)
   - Vectorize loops (compiler pragmas)
   - Cache-friendly memory layouts

### Long-term (Code Rewrite)

7. **Hybrid Python-C++ Architecture**
   - Python for high-level logic and I/O
   - C++/Fortran for performance-critical kernels
   - Use tools like pybind11, f2py
   - Equation specification in symbolic math (SymPy)

8. **Adopt Einstein Toolkit**
   - Industry-standard framework for NR
   - Includes tested modules: BSSN, AMR, analysis
   - Active community support

9. **Auto-generate Code from Equations**
   - Use computer algebra systems (Mathematica, Kranc)
   - Specify equations symbolically
   - Generate optimized C++/Fortran automatically

---

## Core Equations Being Solved

### Einstein Field Equations (EFE)
```
Rμν - ½ gμν R = 8π Tμν
```

### 3+1 Decomposition (ADM Formalism)
- **Metric**: ds² = -α² dt² + γᵢⱼ(dxⁱ + βⁱ dt)(dxʲ + βʲ dt)
- **Extrinsic curvature**: Kᵢⱼ = -½ (∂ₜγᵢⱼ - Dᵢβⱼ - Dⱼβᵢ)
- **Constraints**:
  - Hamiltonian: R + K² - KᵢⱼKⁱʲ = 16π ρ
  - Momentum: Dⱼ(Kⁱʲ - γⁱʲK) = 8π jⁱ

### BSSN Evolution System
17 variables: {φ, K, γ̃ᵢⱼ (6), Ãᵢⱼ (6), Γ̃ⁱ (3), α, βⁱ (3)}

**Key advantage**: Strongly hyperbolic, stable for long-term evolution

### Numerical Methods
- **Spatial discretization**: 2nd/4th/6th/8th order finite differences
- **Time integration**: 4th-order Runge-Kutta (RK4)
- **Boundary conditions**: Sommerfeld outgoing wave
- **Dissipation**: Kreiss-Oliger artificial viscosity
- **AMR**: Berger-Oliger adaptive mesh refinement

---

## Summary

AMSS-NCKU is a **functional but dated** numerical relativity code that successfully simulates binary black hole mergers. Its main strengths are:
- Complete implementation of BSSN/Z4C formulations
- Working AMR with moving grids
- Proven track record (used in Chinese GW research)

Its main weaknesses are:
- Legacy codebase with poor maintainability
- Extensive code duplication
- Manual tensor algebra (error-prone, hard to extend)
- Custom infrastructure instead of standard frameworks

**Recommendation**: For production use, consider migrating to the **Einstein Toolkit**. For educational purposes or specialized modifications, focus on refactoring the RHS modules to use modern tensor libraries.
