# Apical_detachment Project Instructions

This Julia project focuses on simulating and analyzing apical detachment phenomena.

## Project Structure

- `src/`: Main source code directory (contains main module)
- `test/`: Test suite for the project
- `Project.toml`: Julia project manifest with dependencies
- `README.md`: Project documentation
- `.vscode/settings.json`: VS Code Julia configuration

## Getting Started

1. **Install Julia**: Download from [julialang.org](https://julialang.org)
2. **Activate environment**: Open Julia REPL and run:
   ```julia
   cd("path/to/Apical_detachment")
   import Pkg
   Pkg.activate(".")
   ```
3. **Run tests**: In Julia REPL:
   ```julia
   Pkg.test()
   ```

## Development Workflow

- **Code**: Add Julia files to `src/` directory
- **Test**: Add test cases to `test/runtests.jl`
- **Dependencies**: Update `Project.toml` using `Pkg.add("PackageName")`

## IDE Features

The VS Code Julia extension is pre-configured with:
- Julia REPL integration
- Real-time linting and diagnostics
- Code formatting on save
- Debugging support

## Next Steps

1. Add dependencies to `Project.toml` as needed
2. Implement simulation logic in `src/`
3. Add tests to ensure correctness
4. Update this documentation as the project evolves
