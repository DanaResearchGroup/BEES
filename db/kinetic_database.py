#!/usr/bin/env python3

"""
Kinetic Database Module
-----------------------
This module provides access to kinetic parameters from local CSV database with methods for querying enzyme-substrate reactions.

This module knows how to handles:
- Loading kinetic data from CSV files
- Querying by enzyme EC number and substrate
- Handling missing values gracefully
- Logging database access for debugging
"""

import os
import csv
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
from bees.logger import Logger

@dataclass
class KineticData:
    """
    Container for kinetic parameters retrieved from database.
    All concentration units in mM, energy in kJ/mol, temperature in K.
    """
    
    substrate: str
    enzyme: str
    ec_number: str
    cofactors: Optional[str] = None
    products: Optional[str] = None
    amino_acid_sequence: Optional[str] = None  # Enzyme amino acid sequence (uppercase, no spaces)
    km: Optional[float] = None  # mM
    vmax: Optional[float] = None  # mM/s
    kcat: Optional[float] = None  # 1/s
    delta_g: Optional[float] = None  # kJ/mol
    temperature: Optional[float] = None  # K
    ph: Optional[float] = None
    solvent: Optional[str] = None
    inhibitor: Optional[str] = None
    source: str = "database"  # Track data source, not ranking (yet..)
    
    def __repr__(self):
        return (f"KineticData(enzyme={self.enzyme}, substrate={self.substrate}, "
                f"Km={self.km} mM, ΔG={self.delta_g} kJ/mol)")


class KineticDatabase:
    """
    Kinetic parameter database interface.
    
    Loads and queries kinetic data from CSV files. Provides methods to search
    for reactions by enzyme, substrate, and cofactor combinations.
    
    CSV Format Expected:
        index,substrate,cofactors,enzyme,ec_number,amino_acid_sequence,products,temp,km,vmax,kcat,
        delta_g,inhibitor,ph,solvent
    
    Attributes:
        data (List[Dict]): Raw data loaded from CSV
        reactions (List[KineticData]): Parsed kinetic data objects
        logger (Logger): BEES Logger instance for logging database operations
    """
    
    def __init__(self, logger: Logger):
        """
        Initialize empty database.
        
        Args:
            logger (Logger): BEES Logger instance (required).
            
        Raises:
            TypeError: If logger is not a BEES Logger instance.
        """
        if not isinstance(logger, Logger):
            raise TypeError(
                f"logger must be a BEES Logger instance, got {type(logger).__name__}. "
                f"Import from bees.logger and pass a Logger object."
            )
        
        self.data: List[Dict[str, Any]] = []
        self.reactions: List[KineticData] = []
        self.source_file: Optional[str] = None
        self.logger = logger
    
    def load_from_csv(self, csv_path: str) -> int:
        """
        Load kinetic data from CSV file.
        
        Args:
            csv_path (str): Path to CSV file containing kinetic parameters
            
        Returns:
            int: Number of reactions loaded
            
        Raises:
            FileNotFoundError: If CSV file doesn't exist
            ValueError: If CSV format is invalid
        """
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"Kinetic database file not found: {csv_path}")
        
        self.source_file = csv_path
        self.data = []
        self.reactions = []
        
        try:
            with open(csv_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    self.data.append(row)
                    
                    # Parse row into KineticData object
                    try:
                        # Get row index for validation (1-based, since header is row 1)
                        row_index = len(self.reactions) + 2  # +1 for header, +1 for 1-based indexing
                        
                        # Parse and validate parameters
                        km_val = self._parse_float(row.get('km'))
                        vmax_val = self._parse_float(row.get('vmax'))
                        kcat_val = self._parse_float(row.get('kcat'))
                        delta_g_val = self._parse_float(row.get('delta_g'))
                        temp_val = self._parse_float(row.get('temp'))
                        ph_val = self._parse_float(row.get('ph'))
                        
                        # Validate units
                        km_val = self._validate_and_convert_units(km_val, 'km', row_index)
                        vmax_val = self._validate_and_convert_units(vmax_val, 'vmax', row_index)
                        kcat_val = self._validate_and_convert_units(kcat_val, 'kcat', row_index)
                        delta_g_val = self._validate_and_convert_units(delta_g_val, 'delta_g', row_index)
                        temp_val = self._validate_and_convert_units(temp_val, 'temp', row_index)
                        ph_val = self._validate_and_convert_units(ph_val, 'ph', row_index)
                        
                        kinetic_data = KineticData(
                            substrate=row.get('substrate', '').strip(),
                            enzyme=row.get('enzyme', '').strip(),
                            ec_number=row.get('ec_number', '').strip(),
                            cofactors=row.get('cofactors', '').strip() or None,
                            products=row.get('products', '').strip() or None,
                            amino_acid_sequence=row.get('amino_acid_sequence', '').strip() or None,
                            km=km_val,
                            vmax=vmax_val,
                            kcat=kcat_val,
                            delta_g=delta_g_val,
                            temperature=temp_val,
                            ph=ph_val,
                            solvent=row.get('solvent', '').strip() or None,
                            inhibitor=row.get('inhibitor', '').strip() or None,
                            source="database"
                        )
                        self.reactions.append(kinetic_data)
                    except Exception as e:
                        self.logger.warning(f"Failed to parse row in {csv_path}: {row}. Error: {e}")
                        continue
            
            self.logger.info(f"Loaded {len(self.reactions)} kinetic reactions from {csv_path}")
            return len(self.reactions)
            
        except Exception as e:
            raise ValueError(f"Error reading CSV file {csv_path}: {e}")
    
    def _parse_float(self, value: Any) -> Optional[float]:
        """
        Safely parse float from string, handling None, empty, and 'None' strings.
        
        Args:
            value: Value to parse
            
        Returns:
            float or None
        """
        if value is None or value == '' or str(value).strip().lower() == 'none':
            return None
        try:
            return float(value)
        except (ValueError, TypeError):
            return None
    
    def _validate_and_convert_units(self, value: Optional[float], param_name: str, row_index: int) -> Optional[float]:
        """
        Validate that parameter values are in reasonable ranges for their expected units.
        Logs warnings for values that seem outside expected ranges but still returns them.
        
        Expected units:
        - km: mM (typically 0.001 to 1000 mM)
        - vmax: mM/s (typically 0.001 to 1000 mM/s)
        - kcat: 1/s (typically 0.001 to 10000 1/s)
        - delta_g: kJ/mol (typically -200 to 200 kJ/mol)
        - temperature: K (typically 273 to 373 K for biological systems)
        - ph: unitless (must be 0-14)
        
        Args:
            value: Parameter value to validate
            param_name: Name of parameter (km, vmax, kcat, delta_g, temp, ph)
            row_index: Row index in database for error reporting
            
        Returns:
            float or None: Validated value (or None if invalid)
        """
        if value is None:
            return None
        
        # Define reasonable ranges for each parameter
        ranges = {
            'km': (1e-6, 1e6),  # mM: from nanomolar to molar range
            'vmax': (1e-6, 1e6),  # mM/s
            'kcat': (1e-6, 1e7),  # 1/s
            'delta_g': (-500, 500),  # kJ/mol
            'temp': (200, 400),  # K: reasonable biological range
            'ph': (0, 14)  # pH scale
        }
        
        if param_name not in ranges:
            return value  # Unknown parameter, return as-is
        
        min_val, max_val = ranges[param_name]
        
        # Check if value is within reasonable range
        if value < min_val or value > max_val:
            self.logger.warning(
                f"Row {row_index}: {param_name} value {value} may be outside expected range "
                f"[{min_val}, {max_val}] for its units. Please verify units are correct."
            )
        
        # Additional specific checks
        if param_name == 'ph' and (value < 0 or value > 14):
            self.logger.warning(f"Row {row_index}: pH value {value} is outside valid range [0, 14]")
            return None  # pH outside valid range is invalid
        
        if param_name == 'temp' and value < 0:
            self.logger.warning(f"Row {row_index}: Temperature value {value} K is negative (invalid)")
            return None
        
        return value
    
    def query_by_enzyme_substrate(
        self, 
        ec_number: str, 
        substrate_label: str,
        cofactor: Optional[str] = None,
        strict: bool = False,
        temperature_range: Optional[Tuple[float, float]] = None,
        ph_range: Optional[Tuple[float, float]] = None
    ) -> Optional[KineticData]:
        """
        Query database for kinetic parameters by enzyme and substrate.
        
        Args:
            ec_number (str): Enzyme EC number (e.g., "EC 2.7.1.1")
            substrate_label (str): Substrate name/label
            cofactor (str, optional): Cofactor name if required
            strict (bool): If True, require exact cofactor match. If False, ignore cofactor.
            temperature_range (tuple, optional): (min_temp, max_temp) in K. Only return if database
                                                 temperature falls within this range.
            ph_range (tuple, optional): (min_ph, max_ph). Only return if database pH falls within
                                        this range.
            
        Returns:
            KineticData or None: Matching kinetic data or None if not found or outside ranges
        """
        self.logger.debug(f"Querying database: EC={ec_number}, Substrate={substrate_label}, "
                         f"Cofactor={cofactor}, TempRange={temperature_range}, pHRange={ph_range}")
        
        # Normalize inputs
        ec_number = ec_number.strip()
        substrate_label = substrate_label.strip()
        
        for reaction in self.reactions:
            # Check EC number and substrate match
            if (reaction.ec_number == ec_number and 
                reaction.substrate.lower() == substrate_label.lower()):
                
                # Check temperature range if specified
                if temperature_range and reaction.temperature is not None:
                    temp_min, temp_max = temperature_range
                    if not (temp_min <= reaction.temperature <= temp_max):
                        self.logger.debug(f"Database entry temperature {reaction.temperature} K outside "
                                        f"range [{temp_min}, {temp_max}] K")
                        continue
                
                # Check pH range if specified
                if ph_range and reaction.ph is not None:
                    ph_min, ph_max = ph_range
                    if not (ph_min <= reaction.ph <= ph_max):
                        self.logger.debug(f"Database entry pH {reaction.ph} outside "
                                        f"range [{ph_min}, {ph_max}]")
                        continue
                
                # Handle cofactor matching
                if strict and cofactor:
                    if reaction.cofactors and cofactor.lower() in reaction.cofactors.lower():
                        self.logger.debug(f"Database hit: {reaction}")
                        return reaction
                elif not strict:
                    # Return first match regardless of cofactor
                    self.logger.debug(f"Database hit: {reaction}")
                    return reaction
        
        self.logger.debug(f"No database match found for EC={ec_number}, Substrate={substrate_label}")
        return None
    
    def query_by_enzyme(self, ec_number: str) -> List[KineticData]:
        """
        Query all reactions for a given enzyme EC number.
        
        Args:
            ec_number (str): Enzyme EC number
            
        Returns:
            List[KineticData]: All matching reactions
        """
        ec_number = ec_number.strip()
        matches = [rxn for rxn in self.reactions if rxn.ec_number == ec_number]
        self.logger.debug(f"Found {len(matches)} reactions for enzyme {ec_number}")
        return matches
    
    def query_by_substrate(self, substrate_label: str) -> List[KineticData]:
        """
        Query all reactions involving a substrate.
        
        Args:
            substrate_label (str): Substrate name/label
            
        Returns:
            List[KineticData]: All matching reactions
        """
        substrate_label = substrate_label.strip()
        matches = [rxn for rxn in self.reactions 
                  if rxn.substrate.lower() == substrate_label.lower()]
        self.logger.debug(f"Found {len(matches)} reactions for substrate {substrate_label}")
        return matches
    
    def get_all_reactions(self) -> List[KineticData]:
        """
        Get all reactions in the database.
        
        Returns:
            List[KineticData]: All loaded reactions
        """
        return self.reactions
    
    def summary(self) -> Dict[str, Any]:
        """
        Get database summary statistics.
        
        Returns:
            dict: Summary statistics
        """
        total = len(self.reactions)
        with_km = sum(1 for r in self.reactions if r.km is not None)
        with_vmax = sum(1 for r in self.reactions if r.vmax is not None)
        with_kcat = sum(1 for r in self.reactions if r.kcat is not None)
        with_delta_g = sum(1 for r in self.reactions if r.delta_g is not None)
        
        unique_enzymes = len(set(r.ec_number for r in self.reactions))
        unique_substrates = len(set(r.substrate for r in self.reactions))
        
        return {
            'total_reactions': total,
            'unique_enzymes': unique_enzymes,
            'unique_substrates': unique_substrates,
            'reactions_with_km': with_km,
            'reactions_with_vmax': with_vmax,
            'reactions_with_kcat': with_kcat,
            'reactions_with_delta_g': with_delta_g,
            'source_file': self.source_file
        }

