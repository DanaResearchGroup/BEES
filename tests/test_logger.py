import pytest
import os
import shutil
import datetime
import logging
import io
import sys
import time
from typing import Optional, List, Dict, Tuple

# Mock the common module dependencies for testing purposes
class MockMolecule:
    def to_smiles(self):
        return "mock_smiles"

class MockObject:
    def __init__(self, molecule=None, family=None):
        self.molecule = [MockMolecule()] if molecule is None else molecule
        self.family = family

class MockCommon:
    VERSION = "1.0.0"
    bees_path = "/mock/path/to/bees"

    @staticmethod
    def get_git_commit(path):
        return "mock_commit_hash", "mock_commit_date"

    @staticmethod
    def get_git_branch(path):
        return "main"

    @staticmethod
    def time_lapse(t0):
        # Calculate a mock time lapse
        now = datetime.datetime.now()
        delta = now - t0
        total_seconds = int(delta.total_seconds())
        days, remainder = divmod(total_seconds, 86400)
        hours, remainder = divmod(remainder, 3600)
        minutes, seconds = divmod(remainder, 60)
        return f"{days:02d}:{hours:02d}:{minutes:02d}:{seconds:02d}"

    @staticmethod
    def dict_to_str(d):
        return "\n".join([f"{k}: {v}" for k, v in d.items()])

# Inject MockCommon into sys.modules to simulate import
sys.modules['bees.common'] = MockCommon
import bees.common as common # Now this will import our mock

# Actual Logger class code provided by the user, slightly adapted for direct execution in test context
# and to correctly use the mocked common module.
# In a real scenario, you'd just `from src.logger import Logger, logger as global_bees_logger`
# and ensure `bees.common` is truly mockable or correctly setup in your test environment.

logger = logging.getLogger('BEES')
logger.propagate = False
logging.getLogger('matplotlib.font_manager').disabled = True


class Logger(object):
    _initialized = False

    def __init__(self,
                 project: str,
                 project_directory: str,
                 verbose: Optional[int],
                 t0: datetime.datetime,
                 ):

        if Logger._initialized:
            # If the logger has already been initialized, skip reconfiguration to prevent duplicates.
            return

        self.project = project
        self.project_directory = project_directory
        self.t0 = t0
        self.log_file = os.path.join(self.project_directory, 'bees.log') # This seems like an unused attribute in the current Logger code based on review

        self.console_level = verbose if verbose is not None else logging.INFO

        self.main_log_file_path = os.path.join(self.project_directory, f'{self.project}.log')
        self.error_log_file_path = os.path.join(self.project_directory, f'{self.project}_errors.log')

        self.main_file_level = logging.DEBUG
        self.error_file_level = logging.ERROR

        self._setup_handlers()

        Logger._initialized = True

        self.log_header()
        self.info(f"Logger initialized successfully in: {self.project_directory}")
        self.info(f"Console output level: {logging.getLevelName(self.console_level)}")
        self.info(f"Main log file: {self.main_log_file_path} (level: {logging.getLevelName(self.main_file_level)})")
        self.info(f"Error log file: {self.error_log_file_path} (level: {logging.getLevelName(self.error_file_level)})")
        self.info(f"Project '{self.project}' started at: {self.t0.strftime('%Y-%m-%d %H:%M:%S')}")


    def _setup_handlers(self):
        os.makedirs(self.project_directory, exist_ok=True)

        while logger.handlers:
            logger.removeHandler(logger.handlers[0])

        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(self.console_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        if os.path.isfile(self.main_log_file_path):
            self._backup_log_file(self.main_log_file_path, f'{self.project}_main.old')

        main_file_handler = logging.FileHandler(filename=self.main_log_file_path, mode='w')
        main_file_handler.setLevel(self.main_file_level)
        main_file_handler.setFormatter(formatter)
        logger.addHandler(main_file_handler)

        if os.path.isfile(self.error_log_file_path):
            self._backup_log_file(self.error_log_file_path, f'{self.project}_errors.old')

        error_file_handler = logging.FileHandler(filename=self.error_log_file_path, mode='w')
        error_file_handler.setLevel(self.error_file_level)
        error_file_handler.setFormatter(formatter)
        logger.addHandler(error_file_handler)

        logger.setLevel(min(self.console_level, self.main_file_level, self.error_file_level))

    def _backup_log_file(self, file_path: str, base_name: str):
        log_archive_dir = os.path.join(self.project_directory, 'log_archive')
        os.makedirs(log_archive_dir, exist_ok=True)

        local_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        archived_file_name = f'{base_name}.{local_time}.log'
        archived_file_path = os.path.join(log_archive_dir, archived_file_name)

        try:
            shutil.copy(file_path, archived_file_path)
            os.remove(file_path)
            # Use the global logger instance directly for internal logging messages within the Logger class
            global_bees_logger.info(f"Backed up old log file '{os.path.basename(file_path)}' to: {archived_file_path}")
        except Exception as e:
            global_bees_logger.error(f"Failed to backup or remove log file '{file_path}': {e}")


    def debug(self, message: str):
        logger.debug(message)

    def info(self, message: str):
        logger.info(message)

    def warning(self, message: str):
        logger.warning(message)

    def error(self, message: str):
        logger.error(message)

    def critical(self, message: str):
        logger.critical(message)

    def always(self, message: str):
        logger.info(f"ALWAYS: {message}")


    def log_header(self):
        self.always(f'BEES execution initiated on {time.asctime()}\n\n'
                    f'################################################################\n'
                    f'#                                                              #\n'
                    f'#      Biochemical Engine for Enzymatic modelS (BEES)          #\n'
                    f'#                                                              #\n'
                    f'#                                                              #\n'
                    f'# Version:{common.VERSION}{" " * (10 - len(str(common.VERSION)))}          #\n'
                    f'#                                                              #\n'
                    f'################################################################\n\n')

        head,date = common.get_git_commit( path=common.bees_path)
        branch_name = common.get_git_branch(path=common.bees_path)
        if head != '' and date != '':
            self.always(f'The current git HEAD for BEES is:\n'
                                f'    {head}\n    {date}')
        if branch_name and branch_name != 'main':
            self.always(f'    (running on the {branch_name} branch)\n')
        else:
            self.always('\n')
        self.always(f'Starting project {self.project}')


    def log_max_time_reached(self, max_time: str):
        execution_time = common.time_lapse(self.t0)
        self.always(f'Terminating BEES due to time limit.\n'
                    f'Max time set: {max_time}\n'
                    f'Current run time: {execution_time}\n')

    def log_footer(self):
        execution_time = common.time_lapse(self.t0)
        self.always(f'\n\n\nTotal BEES execution time: {execution_time}')
        self.always(f'BEES execution terminated on {time.asctime()}\n')

    def log_species_to_calculate(self,
                                 species_keys: List[int],
                                 species_dict: Dict[int, dict]):
        if len(species_keys):
            self.info('\n\nSpecies to calculate thermodynamic data for:')
            try:
                max_label_length = max([len(spc_dict['QM label']) for key, spc_dict in species_dict.items() if key in species_keys] + [6])
                max_smiles_length = max([len(spc_dict['object'].molecule[0].to_smiles()) for key, spc_dict in species_dict.items() if key in species_keys] + [6])
            except (KeyError, AttributeError):
                self.warning("Could not determine max label/SMILES length for species logging. Missing keys or attributes in species_dict.")
                max_label_length = 20 # Default
                max_smiles_length = 30 # Default

            space1 = ' ' * (max_label_length - len('Label') + 1)
            space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
            self.info(f'Label{space1} SMILES{space2}   Reason for calculating thermo')
            self.info(f'-----{space1} ------{space2}   -----------------------------')
            for key in species_keys:
                spc_dict = species_dict[key]
                try:
                    smiles = spc_dict['object'].molecule[0].to_smiles()
                except (KeyError, AttributeError):
                    smiles = "N/A"
                space1 = ' ' * (max_label_length - len(spc_dict.get('QM label', 'N/A')) + 1)
                space2 = ' ' * (max_smiles_length - len(smiles) + 1)
                reasons = spc_dict.get('reasons', ['No reason provided'])
                one = '1. ' if len(reasons) > 1 else '   '
                self.info(f"{spc_dict.get('QM label', 'N/A')}{space1} {smiles}{space2} {one}{reasons[0]}")
                for j, reason in enumerate(reasons):
                    if j > 0:
                        self.info(f"{' ' * (max_label_length + max_smiles_length + 4)}{j + 1}. {reasons[j]}")


    def log_reactions_to_calculate(self,
                                   reaction_keys: List[int],
                                   reaction_dict: Dict[int, dict]):
        if len(reaction_keys):
            self.info('\n\nReactions to calculate high-pressure limit rate coefficients for:')
            try:
                max_label_length = max([len(rxn_dict['QM label']) for key, rxn_dict in reaction_dict.items() if key in reaction_keys] + [6])
                max_smiles_length = max([len(rxn_dict['SMILES label']) for key, rxn_dict in reaction_dict.items() if key in reaction_keys] + [6])
            except (KeyError, AttributeError):
                self.warning("Could not determine max label/SMILES length for reactions logging. Missing keys or attributes in reaction_dict.")
                max_label_length = 20
                max_smiles_length = 30
            max_label_length = max(max_label_length, max_smiles_length)
            space1 = ' ' * (max_label_length - len('Label') + 1)
            self.info(f'Label{space1} Reason for calculating rate coefficient')
            self.info(f'-----{space1} ---------------------------------------')
            for key in reaction_keys:
                rxn_dict = reaction_dict[key]
                space1 = ' ' * (max_label_length - len(rxn_dict.get('QM label', 'N/A')) + 1)
                self.info(f"\n{rxn_dict.get('QM label', 'N/A')}{space1} {rxn_dict.get('reasons', 'No reason provided')}")
                self.info(f"{rxn_dict.get('SMILES label', 'N/A')}\n")
                if hasattr(rxn_dict.get('object'), 'family') and rxn_dict['object'].family is not None:
                    label = rxn_dict['object'].family if isinstance(rxn_dict['object'].family, str) \
                        else rxn_dict['object'].family.label
                    self.info(f'RMG family: {label}\n')


    def log_species_summary(self, species_dict: Dict[int, dict]):
        if species_dict:
            self.info('\n\n\nSPECIES SUMMARY')
            converged_keys, unconverged_keys = _get_converged_and_unconverged_keys(species_dict)

            try:
                max_label_length = max([len(spc_dict['QM label']) for spc_dict in species_dict.values()] + [6])
                max_smiles_length = max([len(spc_dict['object'].molecule[0].to_smiles()) for spc_dict in species_dict.values()] + [6])
            except (KeyError, AttributeError):
                self.warning("Could not determine max label/SMILES length for species summary. Missing keys or attributes in species_dict.")
                max_label_length = 20
                max_smiles_length = 30

            if len(converged_keys):
                self.info('\nSpecies for which thermodynamic data was calculated:\n')
                space1 = ' ' * (max_label_length - len('label') + 1)
                space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
                self.info(f'Label{space1} SMILES{space2} Reason for calculating thermo for this species')
                self.info(f'-----{space1} ------{space2} ----------------------------------------------')
                for key in converged_keys:
                    spc_dict = species_dict[key]
                    try:
                        smiles = spc_dict['object'].molecule[0].to_smiles()
                    except (KeyError, AttributeError):
                        smiles = "N/A"
                    space1 = ' ' * (max_label_length - len(spc_dict.get('QM label', 'N/A')) + 1)
                    space2 = ' ' * (max_smiles_length - len(smiles) + 1)
                    self.info(f"{spc_dict.get('QM label', 'N/A')}{space1} {smiles}{space2} {spc_dict.get('reasons', 'No reason provided')}")
            else:
                self.info('\nNo species thermodynamic calculation converged!')

            if len(unconverged_keys):
                self.info('\nSpecies for which thermodynamic data did not converge:')
                space1 = ' ' * (max_label_length - len('label') + 1)
                space2 = ' ' * (max_smiles_length - len('SMILES') + 1)
                self.info(f'        Label{space1} SMILES{space2} Reason for calculating thermo for this species')
                self.info(f'        -----{space1} ------{space2} ----------------------------------------------')
                for key in unconverged_keys:
                    spc_dict = species_dict[key]
                    try:
                        smiles = spc_dict['object'].molecule[0].to_smiles()
                    except (KeyError, AttributeError):
                        smiles = "N/A"
                    space1 = ' ' * (max_label_length - len(spc_dict.get('QM label', 'N/A')) + 1)
                    space2 = ' ' * (max_smiles_length - len(smiles) + 1)
                    self.info(f"(FAILED) {spc_dict.get('QM label', 'N/A')}{space1} "
                              f"{smiles}{space2} {spc_dict.get('reasons', 'No reason provided')}")
            else:
                self.info('\nAll species calculated by BEES successfully converged')

    def log_reactions_summary(self, reactions_dict: Dict[int, dict]):
        if reactions_dict:
            self.info('\n\n\nRATE COEFFICIENTS SUMMARY')
            converged_keys, unconverged_keys = _get_converged_and_unconverged_keys(reactions_dict)

            try:
                max_label_length = max([len(rxn_dict['QM label']) for rxn_dict in reactions_dict.values()] + [6])
            except KeyError:
                self.warning("Could not determine max label length for reactions summary. Missing 'QM label' in reaction_dict.")
                max_label_length = 20

            if len(converged_keys):
                self.info('\nReactions for which rate coefficients were calculated:\n')
                space1 = ' ' * (max_label_length - len('label') + 1)
                self.info(f'Label{space1} Reason for calculating rate coefficient for this reaction')
                self.info(f'-----{space1} ---------------------------------------------------------')
                for key in converged_keys:
                    rxn_dict = reactions_dict[key]
                    space1 = ' ' * (max_label_length - len(rxn_dict.get('QM label', 'N/A')) + 1)
                    self.info(f"{rxn_dict.get('QM label', 'N/A')}{space1} {rxn_dict.get('reasons', 'No reason provided')}")
            else:
                self.info('\nNo reaction rate coefficient calculation converged!')

            if len(unconverged_keys):
                self.info('\nReactions for which rate coefficient calculations were unsuccessful:')
                space1 = ' ' * (max_label_length - len('label') + 1)
                self.info(f'        Label{space1} Reason for calculating rate coefficient for this reaction')
                self.info(f'        -----{space1} ---------------------------------------------------------')
                for key in unconverged_keys:
                    rxn_dict = reactions_dict[key]
                    space1 = ' ' * (max_label_length - len(rxn_dict.get('QM label', 'N/A')) + 1)
                    self.info(f"(FAILED) {rxn_dict.get('QM label', 'N/A')}{space1} {rxn_dict.get('reasons', 'No reason provided')}")
            else:
                self.info('\nAll reaction rate coefficients calculations in this iteration successfully converged.')

    def log_unconverged_species_and_reactions(self,
                                              species_keys: List[int],
                                              species_dict: Dict[int, dict],
                                              reaction_keys: List[int],
                                              reaction_dict: Dict[int, dict],
                                              ):
        if len(species_keys):
            self.info('\nThermodynamic calculations for the following species did NOT converge:')
            try:
                max_label_length = max([len(spc_dict['QM label'])
                                        for key, spc_dict in species_dict.items() if key in species_keys] + [6])
            except KeyError:
                self.warning("Could not determine max label length for species logging. Missing 'QM label' in species_dict.")
                max_label_length = 20
            space1 = ' ' * (max_label_length - len('label') + 1)
            self.info(f'Label{space1} SMILES')
            self.info(f'-----{space1} ------')
            for key in species_keys:
                spc_dict = species_dict[key]
                label = spc_dict.get('QM label', 'N/A')
                space1 = ' ' * (max_label_length - len(label))
                try:
                    smiles = spc_dict['object'].molecule[0].to_smiles()
                except (KeyError, AttributeError):
                    smiles = "N/A"
                self.info(f"{label}{space1} {smiles}")
            self.info('\n')
        elif len(species_dict.keys()):
            self.info('\nAll species thermodynamic calculations in this iteration successfully converged.\n')

        if len(reaction_keys):
            self.info('\nRate coefficient calculations for the following reactions did NOT converge:')
            for key in reaction_keys:
                self.info(reaction_dict[key].get('QM label', 'N/A'))
            self.info('\n')
        elif len(reaction_dict.keys()):
            self.info('\nAll reaction rate coefficients calculations in this iteration successfully converged.\n')


    def log_args(self, schema: dict):
        verbose_map = {10: 'debug', 20: 'info', 30: 'warning', None: 'info'}
        schema_copy = schema.copy()
        schema_copy['verbose'] = verbose_map[schema_copy.get('verbose', 20)]

        self.info(f'\n\nUsing the following arguments:\n\n'
                  f'{common.dict_to_str(schema_copy)}')


def _get_converged_and_unconverged_keys(object_dict: Dict[int, dict]) -> Tuple[List[int], List[int]]:
    converged_keys, unconverged_keys = list(), list()
    for key, sub_dict in object_dict.items():
        if 'converged' in sub_dict.keys() and sub_dict['converged']:
            converged_keys.append(key)
        else:
            unconverged_keys.append(key)
    return converged_keys, unconverged_keys


# End of Logger class and supporting functions


# Global logger instance should be accessible as 'global_bees_logger'
global_bees_logger = logging.getLogger('BEES')

# --- Pytest Fixtures ---

@pytest.fixture
def temp_log_dir(tmp_path):
    """
    Creates a temporary directory for log files and ensures it's clean.
    """
    log_dir = tmp_path / "test_logs"
    log_dir.mkdir()
    yield str(log_dir)
    # Clean up after test
    if log_dir.exists():
        shutil.rmtree(log_dir)

@pytest.fixture(autouse=True)
def reset_logger_state():
    """
    Resets the global logger state before each test to ensure tests are independent.
    This is crucial for the Singleton pattern and handler management.
    """
    Logger._initialized = False # Reset the singleton flag
    while global_bees_logger.handlers: # Clear any existing handlers
        global_bees_logger.removeHandler(global_bees_logger.handlers[0])
    global_bees_logger.setLevel(logging.NOTSET) # Reset level
    yield

@pytest.fixture
def capture_console_output():
    """
    Fixture to capture stdout (console output).
    """
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    yield sys.stdout
    sys.stdout = old_stdout

# --- Tests for Logger Class ---

def test_logger_singleton(temp_log_dir):
    """
    Tests that the Logger class correctly implements the Singleton pattern,
    ensuring only one initialization occurs.
    """
    t0 = datetime.datetime.now()
    logger1 = Logger(project="TestProject", project_directory=temp_log_dir, verbose=20, t0=t0)
    assert Logger._initialized is True

    # Attempt to create another instance
    logger2 = Logger(project="AnotherProject", project_directory="/tmp/another_dir", verbose=10, t0=t0)

    # Verify that the second instance did not re-initialize the logger
    # This means logger2 should actually refer to the same configured logger instance (or not change it)
    # The key is that the global_bees_logger.handlers should remain the same as after logger1 init
    # We can check that the project_directory is still the one from logger1
    assert logger2.project_directory == temp_log_dir
    assert global_bees_logger.name == 'BEES'
    assert len(global_bees_logger.handlers) == 3 # Console, main file, error file

    # Ensure the attributes of the second logger instance reflect the first one's configuration
    # because the __init__ was skipped due to _initialized flag
    assert logger2.project == "TestProject"


def test_logger_initialization_and_basic_info_logging(temp_log_dir, capture_console_output):
    """
    Tests logger initialization, header logging, and basic info messages
    to console and main log file.
    """
    t0 = datetime.datetime.now()
    project_name = "TestProject"
    logger_instance = Logger(project=project_name, project_directory=temp_log_dir, verbose=20, t0=t0) # verbose=INFO

    main_log_path = os.path.join(temp_log_dir, f'{project_name}.log')
    error_log_path = os.path.join(temp_log_dir, f'{project_name}_errors.log')

    assert os.path.exists(main_log_path)
    assert os.path.exists(error_log_path)

    logger_instance.info("This is an info message.")
    logger_instance.debug("This is a debug message.") # Should not appear on console (verbose=20)

    console_output = capture_console_output.getvalue()
    assert "This is an info message." in console_output
    assert "This is a debug message." not in console_output # Debug level is lower than info for console

    with open(main_log_path, 'r') as f:
        main_log_content = f.read()
    
    assert "This is an info message." in main_log_content
    assert "This is a debug message." in main_log_content # Main log file logs DEBUG and up


def test_logger_file_backup(temp_log_dir):
    """
    Tests if existing log files are backed up to log_archive upon re-initialization.
    """
    project_name = "BackupTest"
    main_log_path = os.path.join(temp_log_dir, f'{project_name}.log')
    error_log_path = os.path.join(temp_log_dir, f'{project_name}_errors.log')
    archive_dir = os.path.join(temp_log_dir, 'log_archive')

    # Create dummy old log files
    with open(main_log_path, 'w') as f:
        f.write("Old main log content\n")
    with open(error_log_path, 'w') as f:
        f.write("Old error log content\n")

    # Re-initialize logger, which should trigger backup
    t0 = datetime.datetime.now()
    Logger._initialized = False # Reset flag to allow re-initialization
    logger_instance = Logger(project=project_name, project_directory=temp_log_dir, verbose=20, t0=t0)

    assert os.path.exists(archive_dir)
    
    # Check for archived files with timestamp
    archived_main_files = [f for f in os.listdir(archive_dir) if f.startswith(f'{project_name}_main.old.')]
    archived_error_files = [f for f in os.listdir(archive_dir) if f.startswith(f'{project_name}_errors.old.')]
    
    assert len(archived_main_files) == 1
    assert len(archived_error_files) == 1

    # Verify original files are new (or empty if nothing new logged)
    with open(main_log_path, 'r') as f:
        assert "Old main log content" not in f.read()
    with open(error_log_path, 'r') as f:
        assert "Old error log content" not in f.read()


def test_logger_levels(temp_log_dir, capture_console_output):
    """
    Tests that different logging levels are handled correctly by console and file handlers.
    """
    t0 = datetime.datetime.now()
    project_name = "LevelTest"
    logger_instance = Logger(project=project_name, project_directory=temp_log_dir, verbose=logging.WARNING, t0=t0) # Console level: WARNING

    main_log_path = os.path.join(temp_log_dir, f'{project_name}.log')
    error_log_path = os.path.join(temp_log_dir, f'{project_name}_errors.log')

    logger_instance.debug("Debug message")
    logger_instance.info("Info message")
    logger_instance.warning("Warning message")
    logger_instance.error("Error message")
    logger_instance.critical("Critical message")

    console_output = capture_console_output.getvalue()

    # Console output (WARNING level)
    assert "Debug message" not in console_output
    assert "Info message" not in console_output
    assert "Warning message" in console_output
    assert "Error message" in console_output
    assert "Critical message" in console_output

    with open(main_log_path, 'r') as f:
        main_log_content = f.read()

    # Main log file (DEBUG level)
    assert "Debug message" in main_log_content
    assert "Info message" in main_log_content
    assert "Warning message" in main_log_content
    assert "Error message" in main_log_content
    assert "Critical message" in main_log_content

    with open(error_log_path, 'r') as f:
        error_log_content = f.read()

    # Error log file (ERROR level)
    assert "Debug message" not in error_log_content
    assert "Info message" not in error_log_content
    assert "Warning message" not in error_log_content
    assert "Error message" in error_log_content
    assert "Critical message" in error_log_content


def test_logger_always_method(temp_log_dir, capture_console_output):
    """
    Tests the 'always' method and its output format.
    """
    t0 = datetime.datetime.now()
    project_name = "AlwaysTest"
    logger_instance = Logger(project=project_name, project_directory=temp_log_dir, verbose=logging.INFO, t0=t0)

    logger_instance.always("This message should always appear.")

    console_output = capture_console_output.getvalue()
    assert "ALWAYS: This message should always appear." in console_output

    main_log_path = os.path.join(temp_log_dir, f'{project_name}.log')
    with open(main_log_path, 'r') as f:
        main_log_content = f.read()
    assert "ALWAYS: This message should always appear." in main_log_content


def test_logger_log_header_and_footer(temp_log_dir, capture_console_output):
    """
    Tests log_header and log_footer methods.
    Mocks bees.common dependencies for these methods.
    """
    t0 = datetime.datetime.now()
    project_name = "HeaderFooterTest"
    logger_instance = Logger(project=project_name, project_directory=temp_log_dir, verbose=logging.INFO, t0=t0)

    console_output = capture_console_output.getvalue()
    main_log_path = os.path.join(temp_log_dir, f'{project_name}.log')

    # log_header is called during initialization, so it should already be there
    assert "BEES execution initiated" in console_output
    assert "BEES execution initiated" in open(main_log_path).read()
    assert f"Starting project {project_name}" in console_output
    assert f"Starting project {project_name}" in open(main_log_path).read()


    # Test log_footer explicitly
    logger_instance.log_footer()
    console_output_after_footer = capture_console_output.getvalue()
    assert "Total BEES execution time:" in console_output_after_footer
    assert "BEES execution terminated on" in console_output_after_footer
    assert "Total BEES execution time:" in open(main_log_path).read()
    assert "BEES execution terminated on" in open(main_log_path).read()


def test_logger_log_max_time_reached(temp_log_dir, capture_console_output):
    """
    Tests log_max_time_reached method.
    """
    t0 = datetime.datetime.now()
    project_name = "MaxTimeTest"
    logger_instance = Logger(project=project_name, project_directory=temp_log_dir, verbose=logging.INFO, t0=t0)

    max_time_str = "01:00:00:00"
    logger_instance.log_max_time_reached(max_time=max_time_str)

    console_output = capture_console_output.getvalue()
    assert "Terminating BEES due to time limit." in console_output
    assert f"Max time set: {max_time_str}" in console_output
    assert "Current run time:" in console_output

    main_log_path = os.path.join(temp_log_dir, f'{project_name}.log')
    with open(main_log_path, 'r') as f:
        main_log_content = f.read()
    assert "Terminating BEES due to time limit." in main_log_content
    assert f"Max time set: {max_time_str}" in main_log_content
    assert "Current run time:" in main_log_content


def test_logger_log_args(temp_log_dir, capture_console_output):
    """
    Tests log_args method.
    """
    t0 = datetime.datetime.now()
    project_name = "ArgsTest"
    logger_instance = Logger(project=project_name, project_directory=temp_log_dir, verbose=logging.INFO, t0=t0)

    test_schema = {
        "project": "TestProject",
        "verbose": 20,
        "setting1": "value1",
        "setting2": 123
    }
    logger_instance.log_args(test_schema)

    console_output = capture_console_output.getvalue()
    assert "Using the following arguments:" in console_output
    assert "project: TestProject" in console_output
    assert "verbose: info" in console_output # Check verbose mapping
    assert "setting1: value1" in console_output
    assert "setting2: 123" in console_output

    main_log_path = os.path.join(temp_log_dir, f'{project_name}.log')
    with open(main_log_path, 'r') as f:
        main_log_content = f.read()
    assert "Using the following arguments:" in main_log_content
    assert "project: TestProject" in main_log_content
    assert "verbose: info" in main_log_content
    assert "setting1: value1" in main_log_content
    assert "setting2: 123" in main_log_content


def test_logger_verbose_levels(temp_log_dir, capture_console_output):
    """
    Tests the verbose parameter's effect on console output.
    """
    t0 = datetime.datetime.now()

    # Test with verbose=10 (DEBUG)
    Logger._initialized = False # Reset for new logger instance
    logger_debug_console = Logger(project="VerboseDebug", project_directory=temp_log_dir, verbose=10, t0=t0)
    logger_debug_console.debug("Console Debug Message")
    logger_debug_console.info("Console Info Message")
    console_output_debug = capture_console_output.getvalue()
    assert "Console Debug Message" in console_output_debug
    assert "Console Info Message" in console_output_debug
    capture_console_output.truncate(0) # Clear buffer
    capture_console_output.seek(0)

    # Test with verbose=20 (INFO)
    Logger._initialized = False # Reset for new logger instance
    logger_info_console = Logger(project="VerboseInfo", project_directory=temp_log_dir, verbose=20, t0=t0)
    logger_info_console.debug("Console Debug Message")
    logger_info_console.info("Console Info Message")
    console_output_info = capture_console_output.getvalue()
    assert "Console Debug Message" not in console_output_info
    assert "Console Info Message" in console_output_info
    capture_console_output.truncate(0)
    capture_console_output.seek(0)

    # Test with verbose=30 (WARNING)
    Logger._initialized = False # Reset for new logger instance
    logger_warning_console = Logger(project="VerboseWarn", project_directory=temp_log_dir, verbose=30, t0=t0)
    logger_warning_console.debug("Console Debug Message")
    logger_warning_console.info("Console Info Message")
    logger_warning_console.warning("Console Warning Message")
    console_output_warn = capture_console_output.getvalue()
    assert "Console Debug Message" not in console_output_warn
    assert "Console Info Message" not in console_output_warn
    assert "Console Warning Message" in console_output_warn

# Mock data for log_species_to_calculate, log_reactions_to_calculate, log_species_summary, log_reactions_summary
# These mock classes and dictionaries are for testing specific logging methods.

def create_mock_species_dict(converged=True, has_smiles=True, num_reasons=1):
    mock_smiles = MockObject(molecule=[MockMolecule()]) if has_smiles else MockObject(molecule=[])
    reasons = [f"reason {i+1}" for i in range(num_reasons)] if num_reasons > 0 else []
    return {
        1: {
            'QM label': 'SpeciesA',
            'object': mock_smiles,
            'converged': converged,
            'reasons': reasons
        },
        2: {
            'QM label': 'SpeciesB',
            'object': mock_smiles,
            'converged': not converged, # For testing unconverged
            'reasons': ['Another reason']
        }
    }

def create_mock_reaction_dict(converged=True):
    mock_object_with_family = MockObject(family="mock_family_label")
    return {
        101: {
            'QM label': 'Reaction1',
            'SMILES label': 'A=B>>C#D',
            'object': mock_object_with_family,
            'converged': converged,
            'reasons': 'Reason for reaction 1'
        },
        102: {
            'QM label': 'Reaction2',
            'SMILES label': 'E+F>>G',
            'object': mock_object_with_family,
            'converged': not converged,
            'reasons': 'Reason for reaction 2'
        }
    }


def test_log_species_to_calculate(temp_log_dir, capture_console_output):
    t0 = datetime.datetime.now()
    logger_instance = Logger(project="SpeciesCalcTest", project_directory=temp_log_dir, verbose=20, t0=t0)
    
    species_dict = create_mock_species_dict(converged=True, num_reasons=2)
    species_keys = list(species_dict.keys())

    logger_instance.log_species_to_calculate(species_keys, species_dict)

    console_output = capture_console_output.getvalue()
    assert "Species to calculate thermodynamic data for:" in console_output
    assert "Label" in console_output and "SMILES" in console_output
    assert "SpeciesA" in console_output
    assert "SpeciesB" in console_output
    assert "1. reason 1" in console_output # Primary reason
    assert "2. reason 2" in console_output # Secondary reason


def test_log_reactions_to_calculate(temp_log_dir, capture_console_output):
    t0 = datetime.datetime.now()
    logger_instance = Logger(project="ReactionCalcTest", project_directory=temp_log_dir, verbose=20, t0=t0)
    
    reaction_dict = create_mock_reaction_dict(converged=True)
    reaction_keys = list(reaction_dict.keys())

    logger_instance.log_reactions_to_calculate(reaction_keys, reaction_dict)

    console_output = capture_console_output.getvalue()
    assert "Reactions to calculate high-pressure limit rate coefficients for:" in console_output
    assert "Label" in console_output and "Reason for calculating rate coefficient" in console_output
    assert "Reaction1" in console_output
    assert "A=B>>C#D" in console_output
    assert "RMG family: mock_family_label" in console_output


def test_log_species_summary(temp_log_dir, capture_console_output):
    t0 = datetime.datetime.now()
    logger_instance = Logger(project="SpeciesSummaryTest", project_directory=temp_log_dir, verbose=20, t0=t0)
    
    # Test with both converged and unconverged species
    species_dict = {
        1: {'QM label': 'ConvergedSpc', 'object': MockObject(), 'converged': True, 'reasons': 'Converged reason'},
        2: {'QM label': 'FailedSpc', 'object': MockObject(), 'converged': False, 'reasons': 'Failed reason'}
    }

    logger_instance.log_species_summary(species_dict)

    console_output = capture_console_output.getvalue()
    assert "SPECIES SUMMARY" in console_output
    assert "Species for which thermodynamic data was calculated:" in console_output
    assert "ConvergedSpc" in console_output
    assert "Species for which thermodynamic data did not converge:" in console_output
    assert "(FAILED) FailedSpc" in console_output

def test_log_reactions_summary(temp_log_dir, capture_console_output):
    t0 = datetime.datetime.now()
    logger_instance = Logger(project="ReactionSummaryTest", project_directory=temp_log_dir, verbose=20, t0=t0)
    
    # Test with both converged and unconverged reactions
    reactions_dict = {
        101: {'QM label': 'ConvergedRxn', 'object': MockObject(family='A'), 'converged': True, 'reasons': 'Converged reason'},
        102: {'QM label': 'FailedRxn', 'object': MockObject(family='B'), 'converged': False, 'reasons': 'Failed reason'}
    }

    logger_instance.log_reactions_summary(reactions_dict)

    console_output = capture_console_output.getvalue()
    assert "RATE COEFFICIENTS SUMMARY" in console_output
    assert "Reactions for which rate coefficients were calculated:" in console_output
    assert "ConvergedRxn" in console_output
    assert "Reactions for which rate coefficient calculations were unsuccessful:" in console_output
    assert "(FAILED) FailedRxn" in console_output

def test_log_unconverged_species_and_reactions(temp_log_dir, capture_console_output):
    t0 = datetime.datetime.now()
    logger_instance = Logger(project="UnconvergedTest", project_directory=temp_log_dir, verbose=20, t0=t0)
    
    species_dict = {
        1: {'QM label': 'SpcU', 'object': MockObject(), 'converged': False},
        2: {'QM label': 'SpcC', 'object': MockObject(), 'converged': True}
    }
    reaction_dict = {
        101: {'QM label': 'RxnU', 'object': MockObject(family='X'), 'converged': False},
        102: {'QM label': 'RxnC', 'object': MockObject(family='Y'), 'converged': True}
    }

    unconverged_species_keys = [k for k,v in species_dict.items() if not v['converged']]
    unconverged_reaction_keys = [k for k,v in reaction_dict.items() if not v['converged']]

    logger_instance.log_unconverged_species_and_reactions(
        unconverged_species_keys, species_dict,
        unconverged_reaction_keys, reaction_dict
    )

    console_output = capture_console_output.getvalue()
    assert "Thermodynamic calculations for the following species did NOT converge:" in console_output
    assert "SpcU" in console_output
    assert "SpcC" not in console_output # Should not be in unconverged list
    assert "Rate coefficient calculations for the following reactions did NOT converge:" in console_output
    assert "RxnU" in console_output
    assert "RxnC" not in console_output # Should not be in unconverged list