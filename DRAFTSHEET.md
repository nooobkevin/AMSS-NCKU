# DRAFTSHEET - AMSS-NCKU Refactoring Notes

## Tensor Contraction Examples - How to Simplify with einsum

### Current Implementation (Fortran in bssn_rhs.f90)

**Problem**: Raising tensor indices A^ij = g^ik g^jl A_kl (lines 260-282)

```fortran
! Current: 23 lines, error-prone manual indexing
Rxx = gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)

Ryy = gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)

Rzz = gupxz * gupxz * Axx + gupyz * gupyz * Ayy + gupzz * gupzz * Azz + &
      TWO*(gupxz * gupyz * Axy + gupxz * gupzz * Axz + gupyz * gupzz * Ayz)

Rxy = gupxx * gupxy * Axx + gupxy * gupyy * Ayy + gupxz * gupyz * Azz + &
     (gupxx * gupyy + gupxy * gupxy)* Axy + &
     (gupxx * gupyz + gupxz * gupxy)* Axz + &
     (gupxy * gupyz + gupxz * gupyy)* Ayz

Rxz = gupxx * gupxz * Axx + gupxy * gupyz * Ayy + gupxz * gupzz * Azz + &
     (gupxx * gupyz + gupxy * gupxz)* Axy + &
     (gupxx * gupzz + gupxz * gupxz)* Axz + &
     (gupxy * gupzz + gupxz * gupyz)* Ayz

Ryz = gupxy * gupxz * Axx + gupyy * gupyz * Ayy + gupyz * gupzz * Azz + &
     (gupxy * gupyz + gupyy * gupxz)* Axy + &
     (gupxy * gupzz + gupyz * gupxz)* Axz + &
     (gupyy * gupzz + gupyz * gupyz)* Ayz
```

### Proposed Solution 1: Python/NumPy einsum

```python
import numpy as np

# Define tensors as 3x3 arrays (at each grid point)
g_up = np.array([[gupxx, gupxy, gupxz],
                 [gupxy, gupyy, gupyz],
                 [gupxz, gupyz, gupzz]])

A_down = np.array([[Axx, Axy, Axz],
                   [Axy, Ayy, Ayz],
                   [Axz, Ayz, Azz]])

# One line: A^ij = g^ik g^jl A_kl
A_up = np.einsum('ik,jl,kl->ij', g_up, g_up, A_down)

# Extract components
Rxx, Rxy, Rxz = A_up[0, :]
Rxy, Ryy, Ryz = A_up[1, :]
Rxz, Ryz, Rzz = A_up[2, :]
```

**Reduction**: 23 lines of manual tensor indexing → 1 line einsum + 3 lines extraction = ~80% reduction in tensor operation code. Full context (including setup) reduces from 23 lines to ~12 lines.

### Proposed Solution 2: Modern Fortran with MATMUL

```fortran
! Organize data into matrices
real*8, dimension(3,3) :: gup, Adown, Aup, temp

! Fill matrices (symmetric)
gup(1,:) = [gupxx, gupxy, gupxz]
gup(2,:) = [gupxy, gupyy, gupyz]
gup(3,:) = [gupxz, gupyz, gupzz]

Adown(1,:) = [Axx, Axy, Axz]
Adown(2,:) = [Axy, Ayy, Ayz]
Adown(3,:) = [Axz, Ayz, Azz]

! For symmetric tensors: A^ij = g^ik g^jl A_kl
! Due to symmetry: A_kl = A_lk and g^ik = g^ki
! We can write: A^ij = (g^ik A_kl) g^lj
! First contraction: contract g^ik with A_kl over k index
temp = MATMUL(gup, Adown)        ! temp^i_l = g^ik A_kl (sum over k)
! Second contraction: contract temp^i_l with g^jl over l index
Aup = MATMUL(temp, TRANSPOSE(gup))  ! Aup^ij = temp^i_l g^lj = g^ik A_kl g^lj

! Extract results
Rxx = Aup(1,1); Rxy = Aup(1,2); Rxz = Aup(1,3)
Ryy = Aup(2,2); Ryz = Aup(2,3); Rzz = Aup(3,3)
```

**Advantages**:
- Compiler can optimize MATMUL (use BLAS if available)
- Clear mathematical meaning
- Less error-prone

---

## Christoffel Symbol Computation

### Current Implementation (lines 235-257)

```fortran
! Second kind of connection: Γ^i_jk = ½ g^il (∂_j g_lk + ∂_k g_lj - ∂_l g_jk)
! 27 components × 3-8 terms each = ~100 lines
Gamxxx = HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy) + gupxz*(TWO*gxzx - gxxz))
Gamyxx = HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy) + gupyz*(TWO*gxzx - gxxz))
Gamzxx = HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy) + gupzz*(TWO*gxzx - gxxz))
! ... 24 more lines
```

### Proposed Solution: einsum

```python
# First derivatives of metric: dg[i,j,k] = ∂_k g_ij
dg = np.array([[[gxxx, gxyx, gxzx],  # ∂_k g_1j
                [gxyx, gyyx, gyzx],  # ∂_k g_2j
                [gxzx, gyzx, gzzx]], # ∂_k g_3j
               
               [[gxxy, gxyy, gxzy],  # etc.
                [gxyy, gyyy, gyzy],
                [gxzy, gyzy, gzzy]],
               
               [[gxxz, gxyz, gxzz],
                [gxyz, gyyz, gyzz],
                [gxzz, gyzz, gzzz]]])

# Christoffel: Γ^i_jk = ½ g^il (∂_j g_lk + ∂_k g_lj - ∂_l g_jk)
# Build the symmetric combination of derivatives
# dg[i,j,k] = ∂_k g_ij, so we need to rearrange indices:
#   ∂_j g_lk = dg[l,k,j] → use transpose(0, 2, 1) to swap j↔k
#   ∂_k g_lj = dg[l,j,k] → already in correct order
#   ∂_l g_jk = dg[j,k,l] → use transpose(1, 2, 0) for cyclic permutation (i→j, j→k, k→i)
sym_deriv = dg.transpose(0, 2, 1) + dg - dg.transpose(1, 2, 0)
Gamma = 0.5 * np.einsum('il,ljk->ijk', g_up, sym_deriv)

# Extract specific components
Gamxxx = Gamma[0, 0, 0]  # Γ^x_xx
Gamyxx = Gamma[1, 0, 0]  # Γ^y_xx
# etc.
```

**Reduction**: 27 lines → 5 lines (81% reduction)

---

## Code Generation Strategy

### Option 1: Hybrid Python-Fortran

**Workflow**:
1. Write physics equations in Python using SymPy
2. Auto-generate optimized Fortran kernels
3. Call from existing C++ infrastructure via f2py

**Example**:
```python
# physics_equations.py
import sympy as sp
from sympy.tensor import IndexedBase, Idx

# Define symbolic tensors
i, j, k, l = sp.symbols('i j k l', cls=Idx)
g = IndexedBase('g')
A = IndexedBase('A')
gup = IndexedBase('gup')

# Equation: A^ij = g^ik g^jl A_kl
A_up = sp.Sum(sp.Sum(gup[i,k] * gup[j,l] * A[k,l], (k, 0, 2)), (l, 0, 2))

# Generate Fortran code
from sympy.printing.fortran import fcode
print(fcode(A_up, assign_to='Aup(i,j)', source_format='free'))
```

### Option 2: Use Kranc (Mathematica-based code generator)

Kranc is used by Einstein Toolkit to generate Cactus thorns from symbolic equations.

**Workflow**:
1. Specify BSSN equations in Mathematica
2. Kranc generates optimized C code
3. Integrates with Cactus/Carpet AMR framework

**Advantages**:
- Standard in NR community
- Generates cache-optimized code
- Includes boundary conditions, symmetries

### Option 3: Direct Migration to Einstein Toolkit

**Why Einstein Toolkit?**
- Industry-standard framework for numerical relativity
- Includes:
  - McLachlan (BSSN/Z4c thorns)
  - Carpet (AMR)
  - Simfactory (job management)
  - Analysis tools (AHFinderDirect, Multipole, etc.)
- Active community, extensive documentation
- Used by LIGO/Virgo for waveform catalogs

**Migration Path**:
1. Map AMSS-NCKU parameters to Cactus parameter files
2. Use McLachlan for evolution (already has BSSN/Z4c)
3. Use Carpet for AMR (more robust than custom implementation)
4. Use TwoPunctures thorn for initial data (same method!)
5. Validate against AMSS-NCKU results

---

## Specific Files to Refactor (Priority Order)

### High Priority
1. **bssn_rhs.f90** (1700 lines → ~400 with einsum)
2. **Z4c_rhs.f90** (1800 lines → ~400)
3. **bssnEScalar_rhs.f90** (1700 lines → ~400)

### Medium Priority
4. **bssn_class.C** (time evolution loop - add documentation)
5. **TwoPunctures.C** (already reasonable, just needs comments)
6. **cgh.C/patch.C** (AMR - consider replacing with Carpet)

### Low Priority
7. **GPU code** (bssn_gpu.cu - rewrite or remove)
8. **Parallel code** (Parallel.C - consider MPI-3 features)

---

## Example: Ricci Tensor Computation

### Current (lines 377-443)

```fortran
! Ricci tensor R_ij = ∇_k ∇^k γ_ij + ... (many terms)
! Each component computed separately with ~15 lines
call fdderivs(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
Rxx = gupxx * fxx + gupyy * fyy + gupzz * fzz + &
     (gupxy * fxy + gupxz * fxz + gupyz * fyz) * TWO

call fdderivs(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
Ryy = gupxx * fxx + gupyy * fyy + gupzz * fzz + &
     (gupxy * fxy + gupxz * fxz + gupyz * fyz) * TWO
! ... 4 more components
```

### Proposed

```python
# Compute second derivatives: ddg[i,j,k,l] = ∂_k ∂_l g_ij
ddg = compute_second_derivatives(g, grid)

# Ricci: R_ij = g^kl ∂_k ∂_l g_ij - ... (other terms from Christoffels)
R_dd = np.einsum('kl,ijkl->ij', g_up, ddg) + ... # other terms

# Extract components
Rxx, Rxy, Rxz = R_dd[0, :]
Rxy, Ryy, Ryz = R_dd[1, :]
Rxz, Ryz, Rzz = R_dd[2, :]
```

---

## Testing Strategy

1. **Unit Tests**: Test individual tensor operations
   ```python
   def test_raise_indices():
       # Known analytic solution
       g_up = np.eye(3)  # Flat space
       A_down = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
       A_up = raise_indices(g_up, A_down)
       assert np.allclose(A_up, A_down)  # Should be identical in flat space
   ```

2. **Integration Tests**: Compare full RHS computation
   - Run both old and new code on same initial data
   - Verify RHS values agree to machine precision
   - Check conservation laws (ADM mass, etc.)

3. **Evolution Tests**: Compare full simulations
   - GW150914-like system
   - Evolve for ~100M
   - Compare waveforms (mismatch < 1%)
   - Check constraint violations

---

## Performance Considerations

### Current Performance (estimated)
- **Tensor contractions**: ~30-40% of RHS evaluation time
- **Finite differences**: ~40-50%
- **Ghost zone exchanges**: ~10-20%

### Expected Improvements with einsum
- **CPU**: 1.2-1.5× speedup (better cache usage, BLAS acceleration)
- **GPU**: 2-5× speedup (natural fit for tensor cores)
- **Maintainability**: Massive improvement (fewer bugs, easier extensions)

### Why einsum is fast
- Single function call = single kernel launch (GPU)
- Optimizes loop ordering for cache
- Uses BLAS level-3 operations when possible
- JIT compilation (numba) for custom patterns

---

## Next Steps

1. ✅ **Phase 1-4 Analysis Complete**
2. **Create Prototype**:
   - Extract one function (e.g., raise_indices)
   - Implement in Python with einsum
   - Verify against Fortran version
   - Benchmark performance
3. **Scale Up**:
   - Convert full bssn_rhs to Python
   - Create F2Py interface
   - Run regression tests
4. **Documentation**:
   - Write API documentation
   - Create physics tutorial
   - Document equation-to-code mapping

---

## References

1. **BSSN Formulation**: Shibata & Nakamura (1995), Baumgarte & Shapiro (1998)
2. **TwoPunctures**: Ansorg, Brügmann, Tichy (2004)
3. **Einstein Toolkit**: http://einsteintoolkit.org
4. **Kranc**: http://kranccode.org
5. **NumPy einsum**: https://numpy.org/doc/stable/reference/generated/numpy.einsum.html
6. **GRChombo**: https://www.grchombo.org (modern C++ NR framework)
