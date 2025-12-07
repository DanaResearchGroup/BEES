# BEES Test Projects Summary

*This document provides an overview of all test projects in the `/examples/` folder and how to run them*

***How to Run Any Test***
1. open the terminal
2. run the following commends:
```bash
~/BEES
python BEES.py --input_file examples/<TEST_NAME>/input.yml
```
or
```bash
~/BEES
python -m bees.main --input_file examples/<TEST_NAME>/input.yml
```

---

***Test Projects Overview***

Here You can see documantion for all the input test in the project 

### 1. ECClassTest
**Purpose:** Test all major EC enzyme classes to verify proper reaction template generation for different enzyme types.

**Configuration:**
- **Enzymes:** 5 enzymes
  - Oxidoreductase (EC 1.1.1.1)
  - Transferase (EC 2.7.1.1)
  - Hydrolase (EC 3.1.3.1)
  - Lyase (EC 4.2.1.1)
  - Isomerase (EC 5.1.1.1)
- **Substrates:** GenericSubstrate (benzene, C6H6) *be awere that for now cofactors add automaticly if users didn't provide them*
- **Environment:** pH 7.0, 298.15 K (room temp)
- **Expected Reactions:** 5 

**Run Command:**
```bash
python BEES.py --input_file examples/ECClassTest/input.yml
```

---

### 2. DatabaseMatchTest
**Purpose:** Test BEES's ability to find and use kinetic parameters from the database


**Configuration:**
- **Enzymes:** Hexokinase (EC 1.1.1.1)
- **Substrates:** 
  - Glucose (with SMILES)
  - ATP (with SMILES)
- **Environment:** pH 7.4, 310 K (body temperature)
- **Expected Reactions:** 2

**Run Command:**
```bash
python BEES.py --input_file examples/DatabaseMatchTest/input.yml
```

---

### 3. MissingCofactorTest
**Purpose:** Test automatic cofactor addition - verifies BEES automatically adds required cofactors (like ATP) even when they're not in the input species list.

**Configuration:**
- **Enzymes:** Hexokinase (EC 1.1.1.1 - requires ATP)
- **Substrates:** Glucose only (ATP NOT in species list)
- **Environment:** pH 7.4, 310 K
- **Expected Reactions:** 1 (Glucose + ATP → G6P, with ATP auto-added)

**Run Command:**
```bash
python BEES.py --input_file examples/MissingCofactorTest/input.yml
```

---

### 4. DatabaseMissTest
**Purpose:** Test database lookup failure - verifies BEES behavior when no kinetic parameters are found (fallback to template-only generation).

**Configuration:**
- **Enzymes:** GenericEnzyme (EC 2.7.1.1 - Transferase)
- **Substrates:** GenericSubstrate (benzene - not in database)
- **Environment:** pH 7.0, 298.15 K
- **Expected Reactions:** 1 (without kinetics, template only)

**Run Command:**
```bash
python BEES.py --input_file examples/DatabaseMissTest/input.yml
```

---

### 5. TestProject
**Purpose:** Multi-enzyme, multi-substrate test to verify combinatorial reaction generation.

**Configuration:**
- **Enzymes:** 
  - Hexokinase (EC 1.1.1.1)
  - Glutamate dehydrogenase (EC 1.4.1.3)
- **Substrates:** 
  - Glucose
  - L-Glutamate
- **Environment:** pH 7.4, 310 K
- **Expected Reactions:** 4 (2 enzymes × 2 substrates)

**Run Command:**
```bash
python BEES.py --input_file examples/TestProject/input.yml
```

---

### 6. ComprehensiveDemo
**Purpose:** **COMPREHENSIVE DEMONSTRATION** - Showcases ALL BEES reaction generation capabilities in a single test project. This is the ideal test to demonstrate the full range of BEES features.

**Capabilities Demonstrated:**
- **Database Matching:** Glucose + Hexokinase (EC 1.1.1.1) - finds kinetic parameters (Km, kcat, ΔG)
- **Database Miss:** GenericSubstrate + GenericEnzyme (EC 2.7.1.1) - fallback to template-only generation
- **Automatic Cofactor Addition:** Lactose + Lactase (EC 2.7.1.1) - NAD+ cofactor auto-added when missing
-  **All EC Classes:** Tests all 6 major EC classes (1-6)
- **Multi-Enzyme/Substrate:** Multiple combinations demonstrating combinatorial generation
-  **Product Inference:** Reactions without database products use template-based inference

**Configuration:**
- **Enzymes:** 7 enzymes covering all major EC classes
  - Hexokinase (EC 1.1.1.1) - Oxidoreductase, database match
  - GenericEnzyme (EC 2.7.1.1) - Transferase, database miss
  - Lactase (EC 2.7.1.1) - Transferase, missing cofactor
  - ATPase (EC 3.1.3.1) - Hydrolase
  - Invertase (EC 4.2.1.1) - Lyase
  - Isomerase (EC 5.3.1.1) - Isomerase
  - Ligase (EC 6.1.1.1) - Ligase
- **Substrates:** 5 substrates
  - Glucose (database match case)
  - Lactose (missing cofactor case)
  - GenericSubstrate (database miss case)
  - ATP (cofactor)
  - Pyruvate (for oxidoreductase reactions)
- **Environment:** pH 7.4, 310 K (body temperature)
- **Expected Reactions:** Multiple reactions demonstrating all capabilities

**Run Command:**
```bash
python BEES.py --input_file examples/ComprehensiveDemo/input.yml
```

**What to Look For:**
- Reactions with database kinetics (Km, kcat, ΔG values)
- Reactions without kinetics (template-only, parameters need estimation)
- Automatic cofactor addition warnings in logs
- All EC classes represented with appropriate templates
- Product inference for reactions without database products

---

***Test Categories***

### Complexity Tests
- **Minimal:** MissingCofactorTest, DatabaseMissTest (1 reaction each)
- **EC Class Coverage:** ECClassTest (5 reactions)
- **Multi-Component:** TestProject (4 reactions)
- **Comprehensive:** ComprehensiveDemo (all capabilities in one test)

### Database Tests
- **Database Match with Multiple Substrates:** DatabaseMatchTest (2 substrates, tests kinetic parameter retrieval)
- **Database Match Failure:** DatabaseMissTest (1 reaction with no kinetics found)

### Cofactor Tests
- **Automatic Cofactor Addition:** MissingCofactorTest (tests ATP auto-addition when missing from input)

### Enzyme Class Tests
- **EC Class Variety:** ECClassTest (tests all major EC classes)
- **All EC Classes:** ComprehensiveDemo (tests all 6 major EC classes: 1-6)

---

***Output Files***

Each test generates the following files in its project directory:

1. **`reactions_summary.txt`** - Detailed reaction information with kinetic parameters
2. **`<ProjectName>.log`** - Execution log with timestamps
3. **`<ProjectName>_errors.log`** - Error log (if any)
4. **`input.yml`** - Validated input configuration

**Note:** These output files are generated locally when you run examples and are **NOT** tracked in git (they are gitignored). To see what the output format looks like, check the example files in `examples/example_output/` directory.

---

***Quick Test Commands***

```bash
# Run all tests sequentially
cd /home/omerkfir/BEES

python BEES.py --input_file examples/ECClassTest/input.yml
python BEES.py --input_file examples/DatabaseMatchTest/input.yml
python BEES.py --input_file examples/MissingCofactorTest/input.yml
python BEES.py --input_file examples/DatabaseMissTest/input.yml
python BEES.py --input_file examples/TestProject/input.yml
python BEES.py --input_file examples/ComprehensiveDemo/input.yml
```

---

## Database Information

All tests use the kinetic database located at:
```
/home/omerkfir/BEES/db/db.csv
```

---

**Importent Notes*

1. Right now, BEES automatically adds cofactors based on the EC number. in the future we want change this.
2. Right now, BEES does NOT actually run ODE simulations, does NOT estimate parameters using different method, does NOT use different solvers. In the input files they used as a Placeholders for now.


---

---

**Last Updated:** November 2025

