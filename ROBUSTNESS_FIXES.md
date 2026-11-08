# Robustness Fixes Summary

This document describes three robustness fixes implemented to address issues identified in code review.

## Issue 1: Simulation failures may be reported as success

**Problem**: ODE solver failures in the enlarger were logged but didn't propagate to `BEES.execute()` return value or CLI exit code. A failing simulation could be reported as successful.

**Fix**:
1. Added `ode_failed: bool` field to `EnlargerResult` (defaults to `False`)
2. Set `result.ode_failed = True` whenever `sim_result.success` is `False` in enlarger
3. Modified `BEES.execute()` to check `enlarger_result.ode_failed` and return `success: False` when ODE fails
4. Updated `BEES.py` CLI to exit with code 1 when execution results indicate failure
5. Added `database_path` to all result dictionaries for traceability

**Files Modified**:
- `bees/enlarger.py`: Added `ode_failed` field to `EnlargerResult`, set it on ODE failures
- `bees/main.py`: Check `ode_failed` and return early failure result, record `database_path` in all results
- `BEES.py`: Exit with code 1 on failure, improved error reporting

**Backward Compatibility**: Existing code that doesn't check `ode_failed` will still work, but now has access to this explicit failure flag.

---

## Issue 2: Missing thermodynamics silently changes physics

**Problem**: When equilibrator cannot resolve ΔG (e.g., unknown compounds, unbalanced reactions), the code fell back to treating reactions as irreversible with only an INFO-level log message. This silently changed the physics without explicit warning.

**Fix**:
1. Added `thermo_missing_data_policy` setting (schema):
   - `"warn_assume_irreversible"` (default): Log warning and assume irreversible
   - `"silent_assume_irreversible"`: Assume irreversible, log at debug level
   - `"error"`: Fail immediately when thermodynamics cannot be computed
2. Modified `ThermoEngine` to accept and respect `missing_data_policy` parameter
3. Updated `_fallback_thermo()` to log at warning/debug/error level based on policy
4. Enhanced `_log_thermo_fallback_prevalence()` to:
   - Show policy in use
   - Display detailed reaction stoichiometry for first 5 affected reactions
   - Provide clearer guidance on how to resolve (add SMILES, review stoichiometry, or set policy to "error")

**Files Modified**:
- `bees/schema.py`: Added `thermo_missing_data_policy` to `Settings` class
- `bees/thermodynamics.py`: 
  - Added `missing_data_policy` parameter to `ThermoEngine.__init__()`
  - Modified `_fallback_thermo()` to respect policy
  - Updated `compute_keq()` to pass policy through
- `bees/enlarger.py`:
  - Pass `thermo_missing_data_policy` from settings to `ThermoEngine`
  - Enhanced fallback prevalence logging with policy and detailed reaction info
- `projects/commented/input.yml`: Documented new setting
- `README.md`: Added documentation for the new policy

**Backward Compatibility**: Default policy is `"warn_assume_irreversible"`, which maintains similar behavior but with explicit warnings instead of info logs.

---

## Issue 3: Unknown database silently falls back

**Problem**: When the requested database file doesn't exist, `resolve_db_path()` fell back to `db/db.csv` with only a warning. The resolved path wasn't recorded in results, making it unclear which database was actually used.

**Fix**:
1. Added `missing_policy` field to `Database` schema:
   - `"fallback_with_warning"` (default): Fall back to `db/db.csv` with warning
   - `"error"`: Fail immediately if requested database not found
2. Modified `resolve_db_path()` to accept and enforce `missing_policy`
3. Added `database_path` to all execution result dictionaries for traceability
4. Enhanced error message to guide users on policy options

**Files Modified**:
- `bees/schema.py`: Added `missing_policy` to `Database` class
- `bees/main.py`:
  - Updated `resolve_db_path()` signature to accept `missing_policy` parameter
  - Added error handling for `missing_policy="error"`
  - Read policy from database config and pass to `resolve_db_path()`
  - Added `database_path` field to all result dictionaries (iterative and batch modes)
- `projects/commented/input.yml`: Documented new setting
- `README.md`: Added documentation for the new policy

**Backward Compatibility**: Default policy is `"fallback_with_warning"`, maintaining existing behavior.

---

## Testing

New test suite in `tests/test_robustness_fixes.py` covers:

### Issue 1 Tests
- `test_enlarger_result_tracks_ode_failure`: Verify `ode_failed` flag is set
- `test_main_execute_returns_failure_on_ode_failure`: Verify `execute()` returns `success=False`

### Issue 2 Tests
- `test_thermo_engine_respects_error_policy`: Policy="error" raises exception
- `test_thermo_engine_warn_policy_logs_warning`: Policy="warn" logs warning
- `test_thermo_engine_silent_policy_no_warning`: Policy="silent" uses debug log
- `test_enlarger_passes_policy_to_thermo_engine`: Settings propagate correctly

### Issue 3 Tests
- `test_resolve_db_path_error_policy_raises`: Policy="error" raises on missing DB
- `test_resolve_db_path_fallback_policy_warns`: Policy="fallback" warns and continues
- `test_resolve_db_path_succeeds_when_db_exists`: Both policies work when DB exists
- `test_execute_records_database_path_in_results`: Database path recorded in results

### Integration Tests
- `test_missing_db_with_error_policy_fails_before_enlargement`: Early failure on missing DB

---

## Configuration Examples

### Strict Mode (Fail on Any Missing Data)
```yaml
database:
  name: ecoli
  missing_policy: error  # Fail if ecoli.csv doesn't exist

settings:
  thermo_missing_data_policy: error  # Fail if any ΔG cannot be computed
  end_time: 1000.0
```

### Lenient Mode (Continue with Warnings)
```yaml
database:
  name: ecoli
  missing_policy: fallback_with_warning  # Fall back to db.csv with warning

settings:
  thermo_missing_data_policy: warn_assume_irreversible  # Warn but continue
  end_time: 1000.0
```

### Silent Mode (Production/Automated)
```yaml
database:
  name: db
  missing_policy: fallback_with_warning

settings:
  thermo_missing_data_policy: silent_assume_irreversible  # No warnings, log at debug
  end_time: 1000.0
```

---

## Summary of Changes

| Issue | Config Option | Default | Alternative | Breaking? |
|-------|--------------|---------|-------------|-----------|
| ODE failure propagation | (automatic) | - | - | No |
| Missing thermodynamics | `settings.thermo_missing_data_policy` | `"warn_assume_irreversible"` | `"error"`, `"silent_assume_irreversible"` | No |
| Missing database | `database.missing_policy` | `"fallback_with_warning"` | `"error"` | No |

All changes are backward compatible. Existing YAML files will work unchanged with the new defaults maintaining similar behavior (but with improved error reporting).

---

## Files Changed Summary

**Core Changes**:
- `bees/schema.py`: Added `missing_policy` to `Database`, `thermo_missing_data_policy` to `Settings`
- `bees/enlarger.py`: Added `ode_failed` field, enhanced thermo logging
- `bees/main.py`: Early failure checking, database policy enforcement, path recording
- `bees/thermodynamics.py`: Policy-aware fallback behavior
- `BEES.py`: Exit code propagation

**Documentation**:
- `README.md`: Documented new policies
- `projects/commented/input.yml`: Added examples and comments
- `ROBUSTNESS_FIXES.md`: This document

**Tests**:
- `tests/test_robustness_fixes.py`: Comprehensive test suite for all three fixes
