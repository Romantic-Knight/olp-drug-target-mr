# Molecular dynamics (P2/P3) — protocol and server setup

Status: **starting structures built and geometry-checked locally**; the MD itself needs AmberTools +
GROMACS, which are not installed here (checked: no OpenMM, no GROMACS, no antechamber/sqm/tleap
locally, and WSL has no distribution). This document is the method to run once the GPU server is up.

## Decisions this protocol encodes

| Decision | Value | Why |
|---|---|---|
| GDF11 assembly | **monomer** | The published grid box and every docking result are monomer-based, so MD stays comparable. ⚠️ Cys A73 forms a **symmetry-generated interchain disulfide** (SSBOND 5, symmetry 7555), so the biological unit is a covalent dimer; this must be stated as a limitation, not silently ignored. |
| MD start pose | **active-site-directed** | Blind docking put **0 of 3 probe ligands** in the catalytic site, so these trajectories describe a site-directed search, not the receptor's preferred pockets. Must be declared in the methods. |
| Complexes | 4: `GDF11_GW8510`, `SUOX_Sanguinarine`, `SUOX_Spiperone`, `SUOX_Vorinostat` | Updated from the redo: Vorinostat is **no longer** a top SUOX hit (−6.35), so Sanguinarine (−7.29) and Spiperone (−7.24) are carried as the current best. Vorinostat is kept as the paper's comparison case. |
| Protein force field | **AMBER ff19SB** | CHARMM36m is not parameterised for organometallic cofactors, and the SUOX centre is Mo-molybdopterin. This is the technical reason SUOX must take the AMBER route. |
| Water | **TIP3P**, 0.15 M NaCl | Standard, and comparable to the wider MR literature. |
| Ligand force field | **GAFF2 + AM1-BCC** (`antechamber`/`sqm`) — optionally **OpenFF Sage** | GAFF2 is the safer route here because the same AmberTools install handles the metal site. |
| Metal site | **MCPB.py bonded model** (AmberTools) | The only free, reproducible way to give Mo-MPT force-field parameters. |
| Engine | **GROMACS** for MD, **gmx_MMPBSA** for end-state energies | GPU throughput plus a maintained MM-PBSA path from GROMACS trajectories. |

## Protocol

| Stage | Setting |
|---|---|
| Box | dodecahedron, minimum solute-to-edge 1.0-1.2 nm |
| Minimisation | steepest descent, then limited-memory Broyden (to remove the bad contacts tleap can leave) |
| NVT | 100 ps, **position restraints on the ligand** and on protein heavy atoms, solute and solvent **coupled to separate thermostats** (avoids the hot-solvent/cold-solute artefact) |
| NPT | 100 ps, restraints released, 1 bar |
| Production | **3 independent replicates per complex**, different velocity seeds; **100 ns minimum, 200 ns preferred** |
| Integrator | 2 fs, LINCS on H-bonds, PME for long-range electrostatics (4 fs only if hydrogens are re-massaged) |
| Thermostat / barostat | 300 K V-rescale; 1 bar Parrinello-Rahman (C-rescale during equilibration) |
| Convergence | block averaging + autocorrelation time; report **block standard errors, never a single-number error** |

### Why 3 replicates is the floor

A single trajectory cannot distinguish a stable binding mode from one long-lived fluctuation. Three
independent seeds give an honest spread; the docking redo already showed that seed-to-seed spread is
small for the *scoring* (<= 0.64 kcal/mol), and the MD replicates test the same question for the
*dynamics*.

## System size and cost estimate

SUOX complexes are ~2,900 receptor heavy atoms + cofactor + ligand; solvated at 1.0-1.2 nm this is
roughly a **35,000-55,000 atom** system. GROMACS on one modern GPU (RTX 4090 / A100 class) does
about 100-200 ns/day at this size, so

    4 complexes x 3 replicates x 100 ns = 1,200 ns  ->  roughly 6-12 GPU-days

Adding the apo control and any extra ligand roughly doubles it. This is the number to price before
renting.

## Directory layout

```
analysis/md/
  complexes/<name>/{receptor.pdb, ligand.pdb, complex.pdb, contacts.txt}
  complexes/complex_manifest.json     # heavy-atom counts, contacts, clash check
  README.md                           # this file
  pixi.toml                           # AmberTools + GROMACS + analysis environment
  scripts/                            # the run chain, stage by stage
```

`receptor.pdb` deliberately contains **heavy atoms only**. tleap adds hydrogens and picks
protonation states from the force field, which is more reliable than trusting docking hydrogens.

## The Mo-MPT site (read from 2A9A, not assumed)

`analysis/scripts/_probe_mo_coordination.py` lists every atom within 3 A of the Mo, which gives the
first coordination shell directly from the crystal:

| donor | distance to Mo | role |
|---|---|---|
| MTE 1501 S1' | 2.38 A | molybdopterin dithiolene |
| MTE 1501 S2' | 2.43 A | molybdopterin dithiolene |
| CYS 185 SG | 2.41 A | cysteine ligand |
| HOH 7505 O | **1.69 A** | oxo ligand (the short distance identifies it) |
| HOH 7506 O | 2.10 A | hydroxo / water |

Three consequences, each of which corrected an earlier revision of this pipeline:

1. **The two coordinating waters are metal ligands, not solvent, and must be kept.** An earlier
   revision dropped all crystallographic waters (and the Mo itself) to match the docking receptor —
   which would have left MCPB.py with nothing to bond the metal to and a chemically wrong centre.
   `prepare_complex.py` now retains waters within 2.6 A of any metal.
2. **MTE donates two of the five coordinating atoms**, so it cannot be treated as an ordinary
   ligand: Mo and MTE must be parameterised together as one bonded unit.
3. **2A9A contains two sulfates; only res 7503 is catalytic** (5.1 A from the Mo, and its midpoint
   with the Mo is exactly the published grid centre). The surface ion res 7504 (~10.5 A) is a
   crystallisation artefact and is now excluded rather than left in to perturb the electrostatics.

The SUOX MD system is therefore: protein + Mo + MTE + catalytic SO4 + 2 coordination waters +
the docked ligand.

## What has already been verified locally

| Check | Result |
|---|---|
| Complexes assembled in one coordinate frame | 4/4 |
| Steric clashes (heavy-atom pairs < 1.6 A) | **0** in all four |
| Protein contacts within 4 A of the ligand | 10-14 residues per complex |
| SUOX metal site complete | Mo + MTE (24 atoms) + catalytic SO4 (5 atoms) + 2 coordination waters, confirmed present in `receptor.pdb` |
| Surface sulfate excluded | yes - only res 7503 retained |
| Docking protocol proven against a known geometry | SO4 self-docking recovered the crystal site to **0.22 A** |
| Scoring meaningfulness | candidates beat property-matched decoys, AUC 0.747 (P = 1e-4), but 5% of decoys also pass -7 kcal/mol, so no absolute-binding claim is made |

## Known limitations to state in the paper

1. **Mo-MPT is treated with a bonded model** (MCPB.py). It reproduces the site's electrostatics and
   geometry but cannot describe Mo redox chemistry or bond breaking. No chemistry is claimed.
2. **GDF11 is a monomer**, while the crystal implies a covalent dimer.
3. **Site-directed start pose**, not a blind-identified pocket (see the decision table).
4. **MM-PBSA end-state estimates only** — entropic terms are unreliable and the absolute values are
   not interpretable; only same-receptor, same-protocol rankings are.
5. No enhanced sampling (metadynamics etc.). Stated as future work.

## Next steps

1. Stand up the server environment from `pixi.toml` / `scripts/00_setup_server.sh`.
2. Parameterise ligands (`01_param_ligands.sh`) and the Mo-MPT site (`02_param_metal_site.sh`).
3. Build and equilibrate one complex (`03_build_system.sh`, `04_run_md.sh`) and **check the ligand
   RMSD stays bounded before committing to all replicates** — this is the go/no-go for the GPU spend.
4. Only then launch the full 1,200 ns grid.
