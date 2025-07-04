# MSDMS: The Mass Spectrometry Data Management System

## Table of Contents
- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
  - [MSDMS Deployment](#msdms-deployment)
  - [Option A: Upload Spectrum Data via the User Interface](#option-a-upload-spectrum-data-via-the-user-interface)
  - [Option B: Upload spectrum data programmatically](#option-b-upload-spectrum-data-programmatically)
  - [Set up a Python virtual environment](#set-up-a-python-virtual-environment)
- [Accessing the Application](#accessing-the-application)
  - [Deployment context](#deployment-context)
  - [How to connect](#how-to-connect)
- [Explanations](#explanations)
  - [Database schema](#database-schema)
  - [Project structure](#project-structure)
- [Known Issues](#known-issues)
  - [Drag and Drop Fails on Chrome (Ubuntu 22.04)](#drag-and-drop-fails-on-chrome-ubuntu-2204)
---


## Introduction

The Mass Spectrometry Data Management System (MSDMS) is a web application that allows users to create project-specific mass 
spectrometry libraries. By default, MSDMS uses public datasets as a foundation. However, it is possible to add your own data 
in the system to generate customized libraries. MSDMS is composed of 3 entities stored in Docker containers:
- A Python script to retrieve and harmonize data
- A PostgreSQL database to store the data and expose it to the web application
- A JavaScript-based web application to access data, create and export libraries

---

## Prerequisites

Base functionality:
- Docker version 28.0.0 or above

Load data programmatically:
- Python 3
- `pip` version 21 or above
- `virtualenv` version 20 or above

---

## Installation

### MSDMS Deployment

1. Clone the repository:

    ```bash
    git clone https://github.com/lsmsgeneva/MSDMS.git 
    ```

2. Define Database Credentials

    A `.env` file must be created locally (on the same level as the file `docker-compose.yml`) to define the database credentials 
    used by the PostgreSQL container. Required variables are:
    - `POSTGRES_USER`: Database username
    - `POSTGRES_PASSWORD`: Database password
    - `POSTGRES_DB`: Database name
    - `POSTGRES_HOST`: Database host name

    These credentials will be automatically passed to the container on deployment.

3. Build and launch the application:

    In the MSDMS folder, run:
    ```bash
    docker compose up -d
    ```

4. Load Data

    A prepackaged MassBank library is provided in the `data/` folder as a compressed archive (`.tar.gz` format). To load the library, 
    first decompress it:
    ```bash
    tar -xzf data/HMDB.tar.gz -C data/
    tar -xzf data/MSBNK.tar.gz -C data/
    ```
    This will extract the data files into the correct `data/` folder structure required by the system.
    
    If you plan to upload a large dataset (more than 20,000 entries), we recommend using the **programmatic upload method** for 
    better performance and reliability (see [Option B](#option-b-upload-spectrum-data-programmatically)).
    
    > 💡 You can also upload data through the MSDMS web interface after deploying the application (see [Option A](#option-a-upload-spectrum-data-via-the-user-interface)).
    
    **Important:** Keep a backup of your data. Deleting entries through MSDMS will also remove the corresponding files from the 
    `MSDMS_data` folder (used as a Docker volume).

### Option A: Upload Spectrum Data via the User Interface

In MSDMS, you can upload a complete library of spectra directly from the UI by navigating to `Menu > Upload New Data`.

This interface performs the following checks:
1. Validates the format of your MassBank files
2. Verifies that the entries do not already exist in the database
3. Checks for metabolite name conflicts

Preview entries before importing, then click `Write data` to load entries into the database and store the MassBank text files 
in the application container.

> ⚠️ For large datasets (20,000+ entries), we recommend the programmatic upload method for better speed and stability.

### Option B: Upload Spectrum Data Programmatically

1. [Set up a Python virtual environment](#set-up-a-python-virtual-environment)

2. Load your MassBank data library into the database:

    ```bash
    python3 msdms_data_scripts/load_data_in_db.py newData data/
    ```
    You can run this command at any time to add new entries.

3. Update specific database columns efficiently without reloading the entire dataset. For example, to update collision_energy_voltage values:

    ```bash
    python3 msdms_data_scripts/load_data_in_db.py update collision_energy_voltage COLLISION_ENERGY data/
    ```
   - `COLLISION_ENERGY` is the field name in the MassBank file
   - `collision_energy_voltage` is the database column to update
   - `data/` is the folder containing your MassBank files

4. **Optional:** Use the first-level folder names in your data directory as source identifiers:

    ```bash
    python3 msdms_data_scripts/load_data_in_db.py newData data/ --source
    ```

5. **Optional:** Control parallel file loading by specifying the number of workers (default is 8):

    ```bash
    python3 msdms_data_scripts/load_data_in_db.py newData data/ --workers 12
    ```
---

### Set up a Python virtual environment

#### First set-up

1. Install `virtualenv` if it is not already installed:
    ```bash
    pip install virtualenv
    ```
2. Create a virtual environment in your folder of interest:
    ```bash
    virtualenv msdms-env
    ```
3. Activate the virtual environment:
    ```bash
    source msdms-env/bin/activate
    ```
4. Install the required Python libraries using the `requirements.txt` file:
    ```bash
    pip install -r msdms_data_scripts/requirements.txt
    ```

#### Deactivate the virtual environment
```bash
deactivate
```

#### Re-activate virtual environment every time you need to load data programmatically into the database
```bash
source msdms-env/bin/activate
```

---

## Accessing the Application

Once MSDMS is [deployed](#msdms-deployment) using Docker, it runs as a local web service on the host machine.

### Deployment context

MSDMS is designed to be deployed on a local server. It is containerized using Docker and intended to be accessed by other machines 
on the same local network.

### How to connect

You can access the MSDMS web interface by opening a web browser on any machine connected to the same local network and entering 
the IP address of the server in the address bar.

> 💡 Make sure that port `80` is not blocked by the server’s firewall and that the Docker container is running.

---

## Explanations

### Database schema
The database schema is initialized automatically on the first deployment of the container:
1. A schema template is defined in `psql/schema.template.sql`, using `{{DB_OWNER}}` as a placeholder for the database user
2. The script `init-db/00-generate-schema.sh` is executed at startup by PostgreSQL (as part of Docker’s default behavior for scripts inside `/docker-entrypoint-initdb.d`)
3. This script injects the actual credentials from the `.env` file into the template, generating `/tmp/01-init-schema.sql` inside the container
4. The generated SQL is then executed using `psql` to create the initial database schema

This approach ensures that schema setup is both dynamic and reproducible, while avoiding accidental execution of the raw template

### Project structure
The application is split in four containers:
- `nginx`: Server services
- `postgredb`: Contains the PostgreSQL database
- `backend`: Manages the interaction of the frontend with the PostgreSQL database and the data files. Based on `Express` framework
- `frontend`: User interface based on `React` framework

---

## Known Issues

### Drag and Drop Fails on Chrome (Ubuntu 22.04)

There is a known issue with the drag-and-drop functionality on **Chrome running on Ubuntu 22.04**. When attempting to upload files or folders via drag and drop, the behavior is inconsistent and may silently fail or ignore some dropped items. This appears to be a platform-specific limitation of Chrome's file handling on certain Linux environments.

✅ **Workaround**: Use the **Browse** button to upload the folder that contains your MassBank files. This method works reliably across all major browsers and platforms.

🛠️ This issue does not occur on:
- Firefox (Ubuntu 22.04)
- Chrome (Windows/macOS)
- Microsoft Edge
