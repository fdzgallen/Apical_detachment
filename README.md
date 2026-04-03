# Apical_detachment

A Julia project for simulating and analyzing apical detachment phenomena.

## Installation

To use this package, first install Julia from [julialang.org](https://julialang.org). Then, activate the project environment:

```julia
cd("path/to/Apical_detachment")
import Pkg
Pkg.activate(".")
```

## Usage

The main simulation script is `src/1D_Apical_detachment.jl`.

- Edit `src/1D_Apical_detachment.jl` to change mechanical and chemical parameters (e.g., `χ`, `η`, `L`, `k`, `σₐ₀`, `koff`, `kon`, `k0`, `eb`, etc).
- Set initial simulation options at the bottom of the file (partition, Δt, T, output folders, etc).
- Run the simulation from your terminal with Julia:

```bash
cd /home/andreugallen/simulations/Apical_detachment
julia --project=. src/1D_Apical_detachment.jl
```

The script writes output VTU and PNG files under `./VTU/` and `./PNG/` by default.

## Development

To run tests:

```julia
Pkg.test()
```

## License

MIT
