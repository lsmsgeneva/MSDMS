import psycopg2 as psql
from psycopg2 import sql
import re
import os
from tqdm import tqdm
import sys
from rdkit import Chem
from dotenv import load_dotenv
import shutil


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
    dict_spec_info["SPLASH"] = None if len(splash_info)==0 else splash_info[0].split(": ")[1].strip()

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
            cur2.execute(insertMetaboQuery,insertMetaboVariables)
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
    If not, inserts the new spectrum entry linked to a specific metabolite.

    Args:
    -----
    spectrumInfo (dict): A dictionary containing spectral metadata. Expected keys include:
        - ACCESSION, SPLASH, INSTRUMENT_TYPE, DATE, ION_MODE,
          COLLISION_ENERGY, MS_TYPE, SOURCE, FRAGMENTATION_MODE, RETENTION_TIME, FILE_PATH

    metaboliteId (int): The ID of the metabolite in the 'metabolites' table
                        to which this spectrum is related.

    sourceFilePath (str): Full path to the input file (e.g., data/MSBNK/AAFC/...)

    dataPath (str): The base folder for the input data (e.g., data/)

    Returns:
    --------
    None
    """

    # Establish DB connection and create cursor
    connection = msdmsConnection()
    cur = connection.cursor()

    try:
        # Check if spectrum already exists based on accession
        cur.execute("SELECT * FROM spectra WHERE accession = %s;", 
                    [spectrumInfo['ACCESSION']])
        rows = cur.fetchall()
        cur.close()
    except:
        print("INVALID SQL SELECT QUERY")
        connection.rollback()
        return

    # If the spectrum does not exist, insert it
    if len(rows) == 0:
        try:
            cur2 = connection.cursor()
            # Prepare INSERT query and parameters for new spectrum entry
            insertSpectrumQuery = "INSERT INTO spectra( metabolite_id, splash, instrument_type, \
                                collection_date, ionization_mode, collision_energy_voltage,\
                                ms_type, accession, source, fragmentation_method, retention_time, file_path ) \
                                VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"
            insertSpectrumVariables = [
                metaboliteId, spectrumInfo['SPLASH'], spectrumInfo['INSTRUMENT_TYPE'],
                spectrumInfo['DATE'], spectrumInfo['ION_MODE'], spectrumInfo['COLLISION_ENERGY'],
                spectrumInfo['MS_TYPE'], spectrumInfo['ACCESSION'], spectrumInfo['SOURCE'],
                spectrumInfo['FRAGMENTATION_MODE'], spectrumInfo['RETENTION_TIME'], spectrumInfo['FILE_PATH']
            ]

            cur2.execute(insertSpectrumQuery, insertSpectrumVariables)
            connection.commit()

            # --- After successful DB insert, move the file ---
            relative_path = os.path.relpath(sourceFilePath, dataPath)
            destination_file_path = os.path.join('MSDMS_data/massbankFiles', relative_path)

            os.makedirs(os.path.dirname(destination_file_path), exist_ok=True)
            shutil.move(sourceFilePath, destination_file_path)

        except Exception as e:
            print(f"Error inserting spectrum {spectrumInfo['ACCESSION']} or moving file:", e)
            connection.rollback()
        finally:
            connection.close()

def alter_spectra_table(db_col_to_add, ms_field_name, dataPath, db_col_filtering=None, value_of_db_col_fitering=None, type_of_filtering="="):
    """
    Adds new information to a specified column in the 'spectra' table of the msdms database,
    based on values parsed from MassBank text files.

    Args:
    -----
    db_col_to_add (str):
        Name of the column in the 'spectra' table to be updated.

    ms_field_name (str):
        The name of the field in the MassBank file from which to extract the new value
        (e.g., 'AC$MASS_SPECTROMETRY: FRAGMENTATION_METHOD').

    data_path : str
        Path to the root directory containing MassBank .txt files and optional blacklist files.

    db_col_filtering (str, optional):
        Name of a column used to filter which rows to update in the 'spectra' table.

    value_of_db_col_fitering (str, optional):
        Value used in the filtering condition for db_col_filtering.

    type_of_filtering (str, optional):
        Type of comparison to apply when filtering the database column.
        Defaults to "=". Supported operators: '=', '!=', '>', '<', '>=', '<='.

    Example:
    --------
    alter_spectra_table(
        db_col_to_add="fragmentation_method",
        ms_field_name="AC$MASS_SPECTROMETRY: FRAGMENTATION_METHOD",
        dataPath="/home/lsms/Documents/massbank/data/massbankFiles",
        db_col_filtering="ms_type",
        value_of_db_col_fitering="MS",
        type_of_filtering="!="
    )
    """

    # Step 1: If filtering is requested, fetch accessions from the DB
    if db_col_filtering is not None and value_of_db_col_fitering is not None:
        connection = msdmsConnection()
        with connection.cursor() as cur:
            try:
                # Construct the SQL query dynamically with optional filtering
                fetchQuery = sql.SQL(
                    "SELECT accession FROM spectra WHERE {filtering_col} {filtering_type} %s and {col_to_add} is null;"
                ).format(
                    filtering_col=sql.Identifier(db_col_filtering),
                    filtering_type=sql.SQL(type_of_filtering),
                    col_to_add=sql.Identifier(db_col_to_add)
                )
                cur.execute(fetchQuery, (value_of_db_col_fitering,))
                rows = cur.fetchall()
                cur.close()

                # Build the list of filenames to process
                fileList = [elem[0]+".txt" for elem in rows]
            except:
                connection.close()
                raise("INVALID SQL QUERY")    
    else:
        # If no filtering is specified, process all files in the data path
        fileList = os.listdir(dataPath)

    # Step 2: Iterate through files and update the database
    with connection.cursor() as cur2:
        for f in tqdm(fileList):
            fileRows = read_file(os.path.join(dataPath, f))

            # Extract the desired metadata field from the MassBank file
            spectrum_info = [s for s in fileRows if ms_field_name in s]

            if len(spectrum_info) > 0:
                # Get accession ID from file
                accession_filter = [s for s in fileRows if "ACCESSION: " in s][0].strip().split(": ")[1].split("-")[-1]
                # Extract the value to insert (last item on the matched line)
                new_data = spectrum_info[0].strip().split(" ")[-1]

                try:
                    # Dynamically update the specified column in the spectra table
                    update_query = sql.SQL(
                        "UPDATE spectra SET {col_to_add} = %s WHERE accession = %s;"
                    ).format(col_to_add=sql.Identifier(db_col_to_add))

                    cur2.execute(update_query, (new_data, accession_filter))
                except:
                    connection.close()
                    raise("INVALID SQL QUERY")

    # Step 3: Commit changes and close connection
    connection.commit()
    connection.close()
