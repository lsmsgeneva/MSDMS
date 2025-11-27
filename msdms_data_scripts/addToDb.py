import os
import psycopg2 as psql
import re
import shutil
import sys
from dotenv import load_dotenv
from psycopg2 import sql
from psycopg2.extras import execute_values
from rdkit import Chem
from tqdm import tqdm


def read_file(file_path):
    """
    Reads the content of a text file and returns its lines as a list of strings.

    Parameters:
    -----------
    file_path : str
        The absolute or relative path to the text file to be read.

    Returns:
    --------
    list of str
        A list where each element is a line from the file, including newline characters (`\n`).

    Example:
    --------
    >>> lines = read_file("/home/massbankFiles/JP001903.txt")
    >>> print(lines[0])
    'ACCESSION: JP001903\\n'
    """

    try:
        # Open the file in read mode and load all lines into a list
        with open(file_path, "r") as file:
            lines = file.readlines()
        return lines
    except Exception as e:
        # Raise a custom error if the file can't be opened or read
        raise Exception(f"An error occurred while reading the file: {e}")


def parse_compound_info(mb_file):
    """
    Parses compound-level metadata from a MassBank file into a structured dictionary.

    Parameters:
    -----------
    mb_file : list of str
        The content of a MassBank file, provided as a list of text lines.

    Returns:
    --------
    dict
        A dictionary containing parsed compound metadata, including keys like 'ACCESSION',
        'INCHIKEY', 'SMILES', 'EXACT_MASS', and external database identifiers.

    Raises:
    -------
    Exception
        If SMILES parsing fails unexpectedly.

    Example:
    --------
    >>> lines = read_file("/home/massbankFiles/JP001903.txt")
    >>> compound_info = parse_compound_info(lines)
    >>> print(compound_info["ACCESSION"])
    'JP001903'
    """

    # Initialize compound info dictionary with expected keys and default values
    dict_comp_info = {'INCHIKEY': None, 'CAS': None, 'PUBCHEM': None, 'CHEMSPIDER': None, 'CHEBI': None}

    # --- Extract ACCESSION number (last part after the dash) ---
    accession_filter = re.compile(r"^ACCESSION:")
    accession_info = list(filter(accession_filter.match, mb_file))[0].split(": ")[1].split("-")[-1].strip()
    dict_comp_info["ACCESSION"] = accession_info

    # --- Extract all compound-related lines (starting with 'CH$') ---
    compound_filter = re.compile(r"^CH\$")
    compound_info = list(filter(compound_filter.match, mb_file))

    # --- Populate dictionary with parsed compound metadata ---
    for info in compound_info:
        split_info = info.split(": ")
        if split_info[0][3:] == "LINK":
            # CH$LINK entries (e.g., CH$LINK: PUBCHEM 12345)
            split_link = split_info[1].split()
            if split_link[1] not in ["N/A\n", "\n"]:
                dict_comp_info[split_link[0]] = split_link[1]
            else:
                dict_comp_info[split_link[0]] = None
        else:
            # Other CH$ fields (e.g., CH$EXACT_MASS, CH$SMILES)
            attribute_key = split_info[0][3:]
            if attribute_key not in dict_comp_info.keys():
                if split_info[1] not in ["N/A\n", "\n"]:
                    dict_comp_info[attribute_key] = split_info[1].strip()
                else:
                    dict_comp_info[attribute_key] = None

    # --- Normalize SMILES using RDKit (generate canonical version) ---
    if dict_comp_info.get('SMILES'):
        try:
            mol = Chem.MolFromSmiles(dict_comp_info['SMILES'])
            if mol:
                dict_comp_info['SMILES'] = Chem.MolToSmiles(mol, True)
            else:
                print(f"[WARN] Invalid SMILES for accession {dict_comp_info['ACCESSION']}. Setting to empty string.")
                dict_comp_info['SMILES'] = ""
        except Exception as e:
            print(f"[WARN] Failed to parse SMILES for accession {dict_comp_info['ACCESSION']}: {e}")
            dict_comp_info['SMILES'] = ""
    else:
        dict_comp_info['SMILES'] = ""

    # --- Round EXACT_MASS to 4 decimals for consistency ---
    dict_comp_info['EXACT_MASS'] = round(float(dict_comp_info['EXACT_MASS']), 4)

    # --- Ensure COMPOUND_CLASS is always present, even if missing ---
    if 'COMPOUND_CLASS' not in dict_comp_info:
        dict_comp_info['COMPOUND_CLASS'] = None

    return dict_comp_info


def parse_spectra_info(mb_file, file_path, source=""):
    """
    Extracts spectral metadata from a MassBank file and returns it as a structured dictionary.

    Parameters:
    -----------
    mb_file : list of str
        The lines of a MassBank .txt file, typically from `read_file()`.
    file_path : str
        The file path and name of the MassBank .txt file that contains the information to be parsed
    source : str
        Spectrum source. Will be extracted from the accession string if no value is provided. The default value is an
        empty string.

    Returns:
    --------
    dict
        A dictionary containing spectral metadata, including:

    Example:
    --------
    >>> lines = read_file("/data/JP001903.txt")
    >>> spec_info = parse_spectra_info(lines)
    >>> print(spec_info["INSTRUMENT_TYPE"])
    'LC-ESI-QTOF'
    """

    dict_spec_info = {}

    # --- Parse accession ID and infer data source ---
    accession_filter = re.compile(r"^ACCESSION:")
    accession_info = list(filter(accession_filter.match, mb_file))[0].split(": ")[1].strip()
    cleaned_acc = accession_info.split("-")[-1].strip()
    dict_spec_info["ACCESSION"] = cleaned_acc
    if not source:
        acc_upper = accession_info.upper()
        if acc_upper.startswith("HMDB"):
            dict_spec_info["SOURCE"] = "HMDB"
        elif acc_upper.startswith("MSBNK") or acc_upper.startswith("MASSBANK"):
            dict_spec_info["SOURCE"] = "MassBank"
        else:
            dict_spec_info["SOURCE"] = "Custom"
    else:
        dict_spec_info["SOURCE"] = source

    # --- Extract DATE from header (only keep the date, discard time) ---
    date_filter = re.compile(r"^DATE:")
    date_info = list(filter(date_filter.match, mb_file))
    dict_spec_info["DATE"] = date_info[0].split(": ")[1].strip().split()[0]

    # --- SPLASH ID is optional; set to None if not present ---
    splash_filter = re.compile(r"^PK\$SPLASH:")
    splash_info = list(filter(splash_filter.match, mb_file))
    dict_spec_info["SPLASH"] = None if len(splash_info) == 0 else splash_info[0].split(": ")[1].strip()

    # --- Extract instrument and measurement parameters from AC$ lines ---
    spectrum_info_filter = re.compile(r"^AC\$")
    spectrum_info = list(filter(spectrum_info_filter.match, mb_file))

    for info in spectrum_info:
        info_split = info.split(": ")
        if info_split[0] == "AC$MASS_SPECTROMETRY":
            # Example: AC$MASS_SPECTROMETRY: ION_MODE POSITIVE
            mass_spec_info = info_split[1].split()
            dict_spec_info[mass_spec_info[0]] = " ".join(mass_spec_info[1:]).strip()
        elif info_split[0] == "AC$CHROMATOGRAPHY" and info_split[1].startswith("RETENTION_TIME"):
            rt_value = info_split[1].split()[1]
            try:
                dict_spec_info['RETENTION_TIME'] = float(rt_value)
            except ValueError:
                dict_spec_info['RETENTION_TIME'] = None
        else:
            # Example: AC$INSTRUMENT_TYPE: LC-ESI-QTOF
            dict_spec_info[info_split[0][3:]] = info_split[1].strip()

        # --- MS$FOCUSED_ION lines ---
        focused_ion_filter = re.compile(r"^MS\$FOCUSED_ION:")
        focused_ion_info = list(filter(focused_ion_filter.match, mb_file))

        dict_spec_info['CCS'] = None
        dict_spec_info['SV'] = None
        dict_spec_info['COV'] = None
        dict_spec_info['MODIFIER'] = None

        for line in focused_ion_info:
            # --- CCS ---
            if re.match(r"^MS\$FOCUSED_ION:\s*CCS", line):
                # Extract numeric part before 'Å²'
                value_part = re.split(r"CCS[:\s]*", line, maxsplit=1)[-1].strip()
                match = re.search(r"[\d.]+", value_part)
                if match:
                    dict_spec_info['CCS'] = float(match.group(0))

            # --- DMS (SV, COV, Modifier) ---
            elif re.match(r"^MS\$FOCUSED_ION:\s*DMS", line, re.IGNORECASE):
                # Example: MS$FOCUSED_ION: DMS 3000 V,45 V, MeOH
                dict_spec_info['SV'] = None
                dict_spec_info['COV'] = None
                dict_spec_info['MODIFIER'] = None

                value_part = re.split(r"DMS[:\s]*", line, maxsplit=1)[-1].strip()
                parts = [p.strip() for p in value_part.split(",")]  # preserve empties

                numeric_values = []
                for p in parts:
                    match = re.search(r"[-+]?\d*\.\d+|\d+", p)
                    if match:
                        numeric_values.append(float(match.group(0)))
                    else:
                        numeric_values.append(None)

                # Assign SV and COV based on order
                if len(numeric_values) >= 1 and numeric_values[0] is not None:
                    dict_spec_info['SV'] = numeric_values[0]
                if len(numeric_values) >= 2 and numeric_values[1] is not None:
                    dict_spec_info['COV'] = numeric_values[1]

                # Modifier detection
                for p in parts:
                    mod_match = re.search(r"[A-Za-z]{2,5}", p)
                    if mod_match:
                        dict_spec_info['MODIFIER'] = mod_match.group(0).strip()
                        break

    # --- Ensure RETENTION_TIME field exists for schema consistency ---
    if 'RETENTION_TIME' not in dict_spec_info:
        dict_spec_info['RETENTION_TIME'] = None

    # --- Ensure COLLISION_ENERGY field exists for schema consistency ---
    if 'COLLISION_ENERGY' not in dict_spec_info:
        dict_spec_info['COLLISION_ENERGY'] = None

    # --- Ensure FRAGMENTATION_MODE field exists for schema consistency ---
    if 'FRAGMENTATION_MODE' not in dict_spec_info:
        dict_spec_info['FRAGMENTATION_MODE'] = None

    # --- Add the file path and name in the dictionary
    dict_spec_info['FILE_PATH'] = file_path

    return dict_spec_info


def metabolite_query_builder(component_info):
    """
    Dynamically builds a parameterized SQL query to retrieve a metabolite record
    from the 'metabolites' table using a combination of available identifiers.

    Parameters:
    -----------
    component_info : dict
        A dictionary containing metabolite descriptors, typically parsed from MassBank data.
        Required keys include:
            - 'FORMULA', 'NAME'
            - Optional: 'INCHIKEY', 'IUPAC', 'SMILES', 'CAS', 'PUBCHEM', 'CHEBI', 'CHEMSPIDER'

    Returns:
    --------
    tuple
        A tuple of (query_string, query_variables), where:
        - query_string (str): The constructed SQL SELECT statement with placeholders.
        - query_variables (list): A list of values to be passed with the SQL query, in order.

    Example:
    --------
    >>> component_info = {
            'FORMULA': 'C6H12O6',
            'NAME': 'Glucose',
            'INCHIKEY': None,
            'SMILES': 'C(C1C(C(C(C(O1)O)O)O)O)O',
            ...
        }
    >>> query, query_vars = metabolite_query_builder(component_info)
    >>> print(query)
    """

    # Start with mandatory fields: chemical formula and name
    query_vars = [
        component_info['FORMULA'],
        component_info['NAME']
    ]

    query = """SELECT *
               FROM metabolites
               WHERE chemical_formula = %s
                 AND name = %s"""

    # Define optional identifiers and their corresponding DB column names
    columns_to_check = [
        ('INCHIKEY', 'inchikey'), ('IUPAC', 'inchi'),
        ('SMILES', 'smiles'), ('CAS', 'cas_registry_number'),
        ('PUBCHEM', 'pubchem'), ('CHEBI', 'chebi'), ('CHEMSPIDER', 'chemspider')
    ]

    # For each optional identifier, dynamically build SQL filter
    for column, db_column in columns_to_check:
        if component_info[column] is None:
            # Use IS NULL for missing fields to match empty DB entries
            query += f" AND {db_column} IS NULL"
        else:
            # Use parameterized equality for non-null fields
            query += f" AND {db_column} = %s"
            query_vars.append(component_info[column])

    query += ";"

    return query, query_vars


def msdmsConnection():
    """
    Establishes and returns a connection to the 'lsms' PostgreSQL database using credentials from a .env file.

    Returns:
    --------
    conn : psycopg2.extensions.connection
        An active connection object to the local database.

    Raises:
    -------
    Prints an error message and returns None if the connection fails.
    """
    # Load environment variables from .env file
    load_dotenv()

    try:
        conn = psql.connect(
            dbname=os.getenv("POSTGRES_DB"),
            user=os.getenv("POSTGRES_USER"),
            host="localhost",
            password=os.getenv("POSTGRES_PASSWORD"),
            port=5432
        )
        # NOTE: If this script is executed from inside a Docker container,
        #       the host should be changed from 'localhost' to 'postgredb'
        return conn
    except:
        print("Unable to reach the database")


def existInMetaboDb(compoundInfo):
    """
    Checks whether a compound already exists in the 'metabolites' database table.
    If it does not exist, the function inserts the compound and returns its new ID.

    Args:
    -----
    compoundInfo (dict): A dictionary containing compound metadata. Expected keys include:
        - NAME, EXACT_MASS, FORMULA, IUPAC, INCHIKEY, CAS, SMILES,
          PUBCHEM, CHEBI, CHEMSPIDER, ACCESSION

    Returns:
    --------
    int or None:
        The database ID of the compound, either pre-existing or newly inserted.
        Returns None if insertion fails.
    """

    newDataId = None

    # Build SQL query and parameters to check for existing compound
    query, queryVariables = metabolite_query_builder(compoundInfo)

    # Establish DB connection and create cursor
    connection = msdmsConnection()
    cur = connection.cursor()

    try:
        # Execute SELECT query to see if the compound already exists
        cur.execute(query, queryVariables)
        rows = cur.fetchall()
        cur.close()
    except:
        print("INVALID SQL SELECT QUERY")
        connection.rollback()

    # If compound does not exist, insert it into the database
    if len(rows) == 0:
        cur2 = connection.cursor()
        insertMetaboQuery = "INSERT INTO metabolites( name, exact_mass, chemical_formula, \
                            inchi, inchikey, cas_registry_number, smiles, \
                            pubchem, chebi, chemspider ) \
                            VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)\
                            RETURNING id;"
        insertMetaboVariables = [
            compoundInfo['NAME'], compoundInfo['EXACT_MASS'], compoundInfo['FORMULA'],
            compoundInfo['IUPAC'], compoundInfo['INCHIKEY'], compoundInfo['CAS'], compoundInfo['SMILES'],
            compoundInfo['PUBCHEM'], compoundInfo['CHEBI'], compoundInfo['CHEMSPIDER']
        ]

        try:
            # Execute INSERT and retrieve new ID
            cur2.execute(insertMetaboQuery, insertMetaboVariables)
            newDataId = cur2.fetchone()[0]
            connection.commit()
        except:
            print(f"INVALID SQL INSERT QUERY FOR METABO IN {compoundInfo['ACCESSION']}")
            connection.rollback()
        connection.close()
    else:
        # Use existing compound ID
        newDataId = rows[0][0]

    return newDataId


def existInSpectraDb(spectrumInfo, metaboliteId, sourceFilePath, dataPath):
    """
    Checks whether a spectrum record already exists in the 'spectra' table of the database.
    If not, inserts the spectrum along with its associated metadata and validated peak data
    into the corresponding tables.

    During the insertion process, the function verifies the integrity of the MassBank file’s
    peak data. Spectra that lack a PK$PEAK section or contain malformed/non-numeric peak data
    are skipped and not inserted into the database.

    After a successful insertion, the processed MassBank file is moved from its input directory
    to the permanent data storage path under 'MSDMS_data/massbankFiles/'.

    Args:
    -----
    spectrumInfo (dict): A dictionary containing the spectrum metadata. Expected keys include:
        - ACCESSION, SPLASH, INSTRUMENT_TYPE, DATE, ION_MODE,
          COLLISION_ENERGY, MS_TYPE, SOURCE, FRAGMENTATION_MODE,
          RETENTION_TIME, CCS, SV, COV, MODIFIER, FILE_PATH

    metaboliteId (int): The ID of the metabolite in the 'metabolites' table that this spectrum
                        record is associated with.

    sourceFilePath (str): The absolute path to the original MassBank file.

    dataPath (str): The root folder of the dataset (e.g., 'data/'). Used to compute the relative
                    path when moving the processed file to the permanent storage directory.

    Returns:
    --------
    None

    Side Effects:
    -------------
    - Inserts a new record into the 'spectra' table if it does not already exist.
    - Validates and inserts all (m/z, relative intensity) peaks into 'spectrum_peaks'.
    - Skips entries that lack valid or correctly formatted peak data.
    - Commits the transaction upon successful completion.
    - Moves the processed file to:
          MSDMS_data/massbankFiles/<relative_path>
    - Prints diagnostic messages for invalid files, SQL errors, or failed file operations.
    """

    # Establish DB connection and create cursor
    connection = msdmsConnection()
    cur = connection.cursor()

    try:
        cur.execute("SELECT * FROM spectra WHERE accession = %s;", [spectrumInfo['ACCESSION']])
        rows = cur.fetchall()
        cur.close()
    except Exception as e:
        print("INVALID SQL SELECT QUERY:", e)
        connection.rollback()
        return

    if len(rows) == 0:
        try:
            cur2 = connection.cursor()
            insertSpectrumQuery = """
                                  INSERT INTO spectra (metabolite_id, splash, instrument_type,
                                                       collection_date, ionization_mode, collision_energy_voltage,
                                                       ms_type, accession, source, fragmentation_method,
                                                       retention_time, ccs, sv, cov, modifier, file_path)
                                  VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                                          %s) RETURNING spectrum_id; \
                                  """

            insertSpectrumVariables = [
                metaboliteId,
                spectrumInfo.get('SPLASH'),
                spectrumInfo.get('INSTRUMENT_TYPE'),
                spectrumInfo.get('DATE'),
                spectrumInfo.get('ION_MODE'),
                spectrumInfo.get('COLLISION_ENERGY'),
                spectrumInfo.get('MS_TYPE'),
                spectrumInfo.get('ACCESSION'),
                spectrumInfo.get('SOURCE'),
                spectrumInfo.get('FRAGMENTATION_MODE'),
                spectrumInfo.get('RETENTION_TIME'),
                spectrumInfo.get('CCS'),
                spectrumInfo.get('SV'),
                spectrumInfo.get('COV'),
                spectrumInfo.get('MODIFIER'),
                spectrumInfo.get('FILE_PATH')
            ]

            cur2.execute(insertSpectrumQuery, insertSpectrumVariables)
            spectrum_id = cur2.fetchone()[0]

            # Validate peaks before committing
            success = insert_peaks_from_file(connection, spectrum_id, sourceFilePath)
            if not success:
                print(f"Entry {spectrumInfo.get('ACCESSION')} skipped due to invalid or missing peak data.")
                connection.rollback()
                return

            connection.commit()

            relative_path = os.path.relpath(sourceFilePath, dataPath)
            destination_file_path = os.path.join('MSDMS_data/massbankFiles', relative_path)
            os.makedirs(os.path.dirname(destination_file_path), exist_ok=True)
            shutil.move(sourceFilePath, destination_file_path)

        except Exception as e:
            print(f"Error inserting spectrum {spectrumInfo['ACCESSION']} or moving file:", e)
            connection.rollback()
        finally:
            connection.close()


def extract_peaks_from_massbank_file(file_path):
    """
    Parses and validates the PK$PEAK section of a MassBank spectrum file.

    The function extracts triplets of numeric peak values (m/z, intensity, relative intensity)
    from the PK$PEAK section. It ensures that:
      - The section exists in the file.
      - Each peak line has exactly three whitespace-separated columns.
      - All columns contain valid numeric (integer or float) values.

    If any malformed or non-numeric peak line is encountered, the entire entry is considered
    invalid and flagged accordingly.

    Args:
    -----
    file_path (str): Path to the MassBank text file to be parsed.

    Returns:
    --------
    tuple:
        (valid_peaks, has_peak_section, has_valid_peaks)
        where:
        - valid_peaks (list[tuple[float, float]]): A list of (m/z, relative_intensity) pairs.
        - has_peak_section (bool): True if a PK$PEAK section was found in the file.
        - has_valid_peaks (bool): True if all peaks were correctly formatted and numeric.

    Side Effects:
    -------------
    None
    """
    valid_peaks = []
    in_peak_section = False
    has_peak_section = False
    has_valid_peaks = True

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            trimmed = line.strip()

            if trimmed.startswith("PK$PEAK:"):
                in_peak_section = True
                has_peak_section = True
                continue

            if in_peak_section:
                # End of peak section
                if not trimmed or trimmed.startswith("//") or re.match(r"^[A-Z]{2}\$", trimmed):
                    in_peak_section = False
                    continue

                parts = re.split(r"\s+", trimmed)

                # Must have exactly 3 columns
                if len(parts) != 3:
                    has_valid_peaks = False
                    break

                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    rel_int = float(parts[2])
                    valid_peaks.append((mz, rel_int))
                except ValueError:
                    has_valid_peaks = False
                    break

    return valid_peaks, has_peak_section, has_valid_peaks


def insert_peak_batch(connection, spectrum_id, peaks):
    """
    Inserts a batch of validated (m/z, relative intensity) peaks into the 'spectrum_peaks' table.

    This helper function uses PostgreSQL's `execute_values()` for efficient bulk insertion.

    Args:
    -----
    connection (psycopg2.Connection): Active database connection object.
    spectrum_id (int): The ID of the spectrum in the 'spectra' table to which the peaks belong.
    peaks (list[tuple[float, float]]): List of (m/z, relative_intensity) pairs to insert.

    Returns:
    --------
    None

    Side Effects:
    -------------
    - Executes an INSERT operation into the 'spectrum_peaks' table.
    - Does not commit the transaction; commits must be handled by the caller.
    """
    with connection.cursor() as cur:
        query = """
                INSERT INTO spectrum_peaks (spectrum_id, mz_value, rel_int)
                VALUES %s \
                """
        data = [(spectrum_id, mz, rel_int) for mz, rel_int in peaks]
        execute_values(cur, query, data)


def insert_peaks_from_file(connection, spectrum_id, file_path, batch_size=500):
    """
    Validates and inserts all peak data from a MassBank spectrum file into the database.

    This function first calls `extract_peaks_from_massbank_file()` to retrieve and verify
    the peak data. Only spectra that contain a valid PK$PEAK section with properly formatted
    numeric triplets are processed. If the data is missing or malformed, the entry is skipped.

    Peak data is inserted in batches for performance efficiency.

    Args:
    -----
    connection (psycopg2.Connection): Active database connection.
    spectrum_id (int): ID of the corresponding spectrum in the 'spectra' table.
    file_path (str): Path to the MassBank file containing the peak data.
    batch_size (int, optional): Number of peaks to insert per batch (default: 500).

    Returns:
    --------
    bool:
        True  – if all peaks were valid and inserted successfully.
        False – if the file lacked a PK$PEAK section or contained invalid/malformed peak data.

    Side Effects:
    -------------
    - Inserts peak data into the 'spectrum_peaks' table in batch mode.
    - Prints warnings and skips insertion for files with missing or invalid peak data.
    """
    peaks, has_peak_section, has_valid_peaks = extract_peaks_from_massbank_file(file_path)

    if not has_peak_section:
        print(f"Skipping file {file_path}: no PK$PEAK section found.")
        return False
    if not has_valid_peaks:
        print(f"Skipping file {file_path}: invalid or malformed peak data.")
        return False
    if not peaks:
        print(f"Skipping file {file_path}: no valid peaks extracted.")
        return False

    # Insert peaks in batches
    batch = []
    for mz, rel_int in peaks:
        batch.append((mz, rel_int))
        if len(batch) >= batch_size:
            insert_peak_batch(connection, spectrum_id, batch)
            batch = []
    if batch:
        insert_peak_batch(connection, spectrum_id, batch)
    return True
