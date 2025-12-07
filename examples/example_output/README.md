# Example Output Files

This directory contains example output files that demonstrate what BEES generates when you run an example.

## Files in this directory

- **`reactions_summary.example.txt`** - Example reactions summary showing:
  - Reactions with database kinetic parameters (Km, kcat, ΔG, etc.)
  - Reactions without database parameters (template-only)
  - Different EC classes and reaction types
  - Cofactor information
  - Rate law types

## What gets generated when you run an example?

When you run any example (e.g., `python BEES.py --input_file examples/ComprehensiveDemo/input.yml`), BEES will generate:

1. **`<ProjectName>.log`** - Main execution log with timestamps and detailed information
2. **`<ProjectName>_errors.log`** - Error log (if any errors occur)
3. **`reactions_summary.txt`** - Summary of all generated reactions with kinetic parameters
4. **`log_archive/`** - Directory containing archived old log files (if you run the same project multiple times)

## Note

These output files are **NOT** tracked in git (they are gitignored) because they are generated locally when you run examples. This directory (`example_output/`) contains example files to show you what the output format looks like.

## Example Output Location

When you run an example, the output files are generated in the same directory as the input file. For example:
- Input: `examples/ComprehensiveDemo/input.yml`
- Output: `examples/ComprehensiveDemo/ComprehensiveDemo.log`, `examples/ComprehensiveDemo/reactions_summary.txt`, etc.

These output files will be in your local filesystem but will not be pushed to GitHub.


