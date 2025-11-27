# -*- coding: utf-8 -*-
"""
MassBank Data Loader and Updater Script
---------------------------------------

This script processes MassBank-formatted `.txt` files to insert new compound and spectra
data into a database or update existing entries. It supports recursive directory scanning,
blacklist filtering, parallel file processing, and selective field updates in the 'spectra' table.
When using the new data functionality, associated spectrum peak information is also extracted
and stored, but peak data are not affected by update operations.

Functions:
----------
- load_data(data_path, new=True, use_folder_as_source=False, num_workers=8):
    Scans for MassBank `.txt` files in the provided path, skips blacklisted files,
    parses compound and spectra information, and inserts it into the database using
    multiprocessing for efficiency.

- updateSpectraTable(dbCol, fieldName, dataPath, use_folder_as_source=False):
    Recursively scans for `.txt` files in a directory, compares a specified field
    in the parsed spectra against existing DB values, and updates entries if mismatches are found.

Dependencies:
-------------
- Python standard libraries: os, sys, re, glob
- tqdm for progress indication
- multiprocessing for parallel file processing
- Custom module: `addToDb` (provides DB connection, file parsing, and insertion logic)

Usage (via command line):
-------------------------
$ python3 load_data_in_db.py newData <data_directory> [--source] [--workers N]
    → Loads all `.txt` files in the directory (recursively) into the database,
      excluding any blacklisted files. If `--source` is provided, the `source`
      field in the spectra table is set to the name of the first-level folder
      containing each file. Use `--workers` to specify the number of parallel
      processes (default is 8).

$ python3 load_data_in_db.py update <db_column> <field_name> <data_directory> [--source]
    → Updates an existing column in the database if the parsed value in the field
      differs from the current DB entry. If `--source` is provided, the `source`
      field in the spectra table is interpreted as the name of the first-level
      folder containing each file.

Examples:
---------
$ python3 msdms_data_scripts/load_data_in_db.py newData data/
$ python3 msdms_data_scripts/load_data_in_db.py newData data/ --source
$ python3 msdms_data_scripts/load_data_in_db.py newData data/ --workers 12
$ python3 msdms_data_scripts/load_data_in_db.py newData data/ --source --workers 16

$ python3 msdms_data_scripts/load_data_in_db.py update collision_energy_voltage COLLISION_ENERGY data/
$ python3 msdms_data_scripts/load_data_in_db.py update collision_energy_voltage COLLISION_ENERGY data/ --source

Notes:
------
- The blacklist mechanism uses any `legacy.blacklist` files found in subdirectories.
- The `load_data` function uses multiprocessing for faster data insertion.
- Be cautious when running the update function — incorrect field mapping may cause SQL errors.
"""

import glob
import os
import re
import sys
from multiprocessing import Pool
from tqdm import tqdm

import addToDb


def process_single_file(args):
    """
    Process a single MassBank .txt file: parse, filter, insert into the database, and move.

    Parameters:
    -----------
    args : tuple
        A tuple containing:
        - f (str): Path to the .txt file.
        - black_listed (list): List of blacklisted filenames.
        - data_path (str): Root directory path.
        - use_folder_as_source (bool): Whether to use the folder name as source.

    Returns:
    --------
    None
    """
    f, black_listed, data_path, use_folder_as_source = args
    try:
        file_name = os.path.basename(f)
        file_name_filter = re.compile(re.escape(file_name))
        black_listed_element = list(filter(file_name_filter.match, black_listed))

        if len(black_listed_element) == 0:
            file_rows = addToDb.read_file(f)
            comp_info = addToDb.parse_compound_info(file_rows)
            relative_path = os.path.relpath(f, data_path)

            if use_folder_as_source:
                first_layer_folder = relative_path.split(os.sep)[0]
                spec_info = addToDb.parse_spectra_info(file_rows, relative_path, source=first_layer_folder)
            else:
                spec_info = addToDb.parse_spectra_info(file_rows, relative_path)

            compound_id = addToDb.existInMetaboDb(comp_info)
            addToDb.existInSpectraDb(spec_info, compound_id, f, data_path)
    except Exception as e:
        print(f"Error processing file {f}: {e}")


def load_data(data_path, use_folder_as_source=False, num_workers=8):
    """
    Load and process MassBank .txt files in parallel, inserting valid entries into the database.

    Recursively scans a directory for .txt files (excluding blacklisted ones), parses compound and
    spectra information, and adds them to the database. Files are processed concurrently using a
    multiprocessing pool.

    Parameters:
    -----------
    data_path : str
        Root directory containing MassBank .txt files and optional legacy.blacklist files.
    use_folder_as_source : bool, optional (default=False)
        If True, sets the 'source' field to the top-level folder name of the file.
    num_workers : int, optional (default=8)
        Number of parallel processes for file parsing and database insertion.

    Returns:
    --------
    None
    """

    # --- Step 1: Collect all .txt files and blacklist files ---
    txt_files = glob.glob(os.path.join(data_path, '**/*.txt'), recursive=True)
    blacklist = glob.glob(os.path.join(data_path, '**/legacy.blacklist'), recursive=True)

    # --- Step 2: Read and accumulate all blacklisted filenames ---
    black_listed = []
    for bf in blacklist:
        black_listed_files = addToDb.read_file(bf)
        black_listed += black_listed_files

    # Prepare arguments for pool
    pool_args = [(f, black_listed, data_path, use_folder_as_source) for f in txt_files]

    # Process files in parallel
    with Pool(processes=num_workers) as pool:
        list(tqdm(pool.imap_unordered(process_single_file, pool_args), total=len(pool_args)))


def updateSpectraTable(dbCol, fieldName, dataPath, use_folder_as_source=False):
    """
    Updates a specified column in the 'spectra' table based on new data
    parsed from MassBank .txt files located recursively in the given path.

    Parameters:
    -----------
    dbCol : str
        The column name in the database to update (e.g., 'INSTRUMENT_TYPE').

    fieldName : str
        The corresponding key in the parsed spectrum dictionary
        (usually same as dbCol but in uppercase).

    dataPath : str
        Root directory containing the MassBank .txt files (searched recursively).

    use_folder_as_source : bool, optional (default=False)
        If True, passes the top-level folder name as 'source' to the spectrum parser.

    Returns: None
    --------
    """

    # --- Step 1: Recursively collect all .txt files ---
    file_list = glob.glob(os.path.join(dataPath, '**/*.txt'), recursive=True)

    # --- Step 2: Fetch current values from DB for comparison ---
    connection = addToDb.msdmsConnection()
    cur = connection.cursor()
    fetchQuery = addToDb.sql.SQL("SELECT {field}, accession FROM spectra;").format(
        field=addToDb.sql.Identifier(dbCol)
    )
    cur.execute(fetchQuery)
    rows = cur.fetchall()
    connection.close()

    accessions = [row[1] for row in rows]

    # --- Step 3: Check and update mismatched entries ---
    for file_path in tqdm(file_list):
        fileRows = addToDb.read_file(file_path)
        relative_path = os.path.relpath(file_path, dataPath)

        # Get top-level folder if needed
        if use_folder_as_source:
            source = relative_path.split(os.sep)[0]
            specInfo = addToDb.parse_spectra_info(fileRows, relative_path, source=source)
        else:
            specInfo = addToDb.parse_spectra_info(fileRows, relative_path)

        accession = specInfo.get("ACCESSION")
        if accession in accessions:
            db_value = rows[accessions.index(accession)][0]
            new_value = specInfo.get(fieldName)

            if new_value != db_value:
                connection = addToDb.msdmsConnection()
                cur2 = connection.cursor()
                updateQuery = addToDb.sql.SQL(
                    "UPDATE spectra SET {field} = %s WHERE accession = %s;"
                ).format(field=addToDb.sql.Identifier(dbCol))

                try:
                    cur2.execute(updateQuery, [new_value, accession])
                    connection.commit()
                except:
                    print(f"INVALID SQL UPDATE QUERY FOR SPECTRA IN {accession}")
                    connection.rollback()
                connection.close()


args = sys.argv

# Check for --help or -h flags
if '--help' in args or '-h' in args:
    print(__doc__)
    sys.exit()

if len(args) < 2:
    print("Insufficient arguments provided. Use --help for usage information.")
    sys.exit()

if args[1] == "newData":
    use_source = '--source' in args

    # Extract --workers value if provided
    num_workers = 8  # Default value
    if '--workers' in args:
        try:
            num_workers_idx = args.index('--workers') + 1
            num_workers = int(args[num_workers_idx])
        except (IndexError, ValueError):
            print("Invalid or missing value for --workers. Using default of 8.")

    load_data(args[2], use_folder_as_source=use_source, num_workers=num_workers)

elif args[1] == "update":
    if len(args) < 5:
        print("Missing required arguments for 'update' mode. Use --help for usage information.")
        sys.exit()
    use_source = '--source' in args
    updateSpectraTable(args[2], args[3], args[4], use_folder_as_source=use_source)

else:
    print(f"Unrecognized command '{args[1]}'. Use --help for usage information.")
