from eHive import BaseRunnable
import os
from pathlib import Path
import sys
import pymysql
import logging

# Import path setup
script_dir = str(Path(__file__).resolve().parents[3]/ "ensembl-genes"/ "src" / "python" / "ensembl" / "genes" / "info_from_registry")
sys.path.append(script_dir)

from build_anno_commands import build_annotation_commands
from check_if_annotated import check_if_annotated
from create_pipe_reg import create_registry_entry
from mysql_helper import mysql_fetch_data
from assign_clade_based_on_tax import (
    assign_clade,
    assign_clade_info_custom_loading,
)

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
if not logger.handlers:  # Avoid duplicate handlers
    stdout_handler = logging.StreamHandler()
    stdout_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stdout_handler.setFormatter(formatter)
    logger.addHandler(stdout_handler)


class ProcessGCApython(BaseRunnable):
    """Process metadata data for the annotation pipeline."""

    def fetch_input(self):
        """Fetch and prepare input data for the annotation pipeline."""
        try:
            assembly_accession = self.param('assembly_accession')
            logger.info(f"Processing assembly accession: {assembly_accession}")

            production_gca = assembly_accession.replace('.', 'v').replace('_', '').lower()
            logger.debug(f"Generated production GCA: {production_gca}")

            core_dbname = f"{self.param('dbowner')}_{production_gca}_core_{self.param('ensembl_release')}_1"
                 
            # Prepare server info
            server_info = {
                "registry": {
                    "db_host": self.param("registry_db")["-host"],
                    "db_user": "ensro",
                    "db_port": self.param("registry_db")["-port"],
                    "db_password": ""
                },
                "pipeline_db": {
                    "db_host": self.param("pipe_db")["-host"],
                    "db_user": self.param("pipe_db")["-user"],
                    "db_port": self.param("pipe_db")["-port"],
                    "db_password":self.param("pipe_db")["-pass"],
                    "db_name": self.param("pipe_db")["-dbname"],
                },
                "core_db": {
                    "db_host": self.param("core_db")["-host"],
                    "db_user": self.param("core_db")["-user"],
                    "db_port": self.param("core_db")["-port"],
                    "db_password": self.param("core_db")["-pass"],
                    "db_name": core_dbname,
                }
            }

            # Handle custom initialization file if provided
            init_file = self.param("init_file")
            if init_file is not None and Path(init_file).is_file():
                logger.info("Init file detected, checking if annotated.")
                check_if_annotated(assembly_accession, server_info)


            # Add data from registry and create a dictionary
            clade_settings_path = self.param("clade_settings_path")
            info_dict = add_generated_data(server_info, assembly_accession, init_file, clade_settings_path, production_gca)
            if info_dict is None:
                raise ValueError("registry_info is None") 

            logger.info("Added metadata from the registry")
        
        
            info_dict["core_db"] = server_info["core_db"]
            base_output_dir = Path(self.param("base_output_dir"))
            output_path = base_output_dir / assembly_accession
            short_read_dir = output_path / "short_read_fastq"

            use_existing_dir = self.param("use_existing_short_read_dir")

            if use_existing_dir != "No" and os.path.isdir(use_existing_dir):
                short_read_dir = use_existing_dir

            
            repeatmodeler_library = self.param("repeatmodeler_library")
            if repeatmodeler_library == "No" or not os.path.isdir(repeatmodeler_library):
                repeatmodeler_library = ""

            # Add values back to dictionary
            info_dict.update(
                {
                    "ensembl_release": self.param("ensembl_release"),
                    "output_path": output_path / assembly_accession,
                    "genome_files_dir": output_path / "genome_files",
                    "toplevel_genome_file": output_path / f"{info_dict['species_name']}_toplevel.fa",
                    "reheadered_toplevel_genome_file": output_path / f"{info_dict['species_name']}_reheadered_toplevel.fa",
                    "short_read_dir": str(short_read_dir),
                    "long_read_dir": output_path / "long_read_fastq",
                    "gst_dir": output_path / "gst",
                    "rnaseq_summary_file": output_path / f"{info_dict['production_name']}.csv",
                    "rnaseq_summary_file_genus": output_path / f"{info_dict['production_name']}_gen.csv",
                    "long_read_summary_file": output_path / "long_read_fastq" / f"{info_dict['production_name']}_long_read.csv",
                    "repeatmodeler_library": repeatmodeler_library,
                    "current_genebuild": self.param("current_genebuild"),
                    "assembly_accession": assembly_accession,
                    "core_dbname": self.param("core_dbname"),
                    "num_threads": self.param("num_threads"),
                }
            )



            # Create directories
            create_dir(info_dict["output_path"], mode=0o775)
            for dir_path in [
                info_dict["genome_files_dir"],
                info_dict["short_read_dir"],
                info_dict["long_read_dir"],
                info_dict["gst_dir"],
            ]:
                create_dir(dir_path)

            logger.info("Directories created")

            # Registry and annotation commands
            adaptors = {
                "core_string": {
                    "host": server_info["core_db"]["db_host"],
                    "port": server_info["core_db"]["db_port"],
                    "dbname": self.param("core_dbname"),
                    "user": self.param("user"),
                    "pass": self.param("password"),
                    "species": info_dict["production_name"],
                    "group": "core",
                },
            }
        
            registry_path = create_registry_entry(server_info, adaptors, base_output_dir, self.param("registry_file"))

            info_dict["registry_file"] = str(registry_path)
            build_annotation_commands(adaptors, info_dict)
            logger.debug("Anno commands created successfully")

            # Store output for write_output
            self.output_params = info_dict
            logger.info("Input fetch completed successfully")

        except Exception as e:
            logger.error(f"Error in fetch_input: {e}")
            raise

    def run(self):
        """Run method - currently empty as processing is done in fetch_input."""
        pass


    def write_output(self):
        self.dataflow(self.output_params, 1)
        


from pathlib import Path

def get_metadata_from_registry(server_info, assembly_accession, init_file):
    """
    Retrieves registry metadata for a given genome assembly accession.

    Supports two modes:
    1. Use custom registry file (`init_file`)
    2. Query MySQL registry database
    """
    # Use custom registry mode
    if init_file is not None:
        p = Path(init_file)
        if p.is_file():
            registry_info = custom_loading(p)
            registry_info = assign_clade_info_custom_loading(registry_info)
            return registry_info
        else:
            raise FileNotFoundError(f"INI file not found: {init_file}")

    # Fallback to MySQL registry
    try:
        assembly_accessions = (
            [assembly_accession] if isinstance(assembly_accession, str) else assembly_accession
        )
        placeholders = ",".join(["%s"] * len(assembly_accessions))

        registry_query = f"""
            SELECT 
                s.species_taxon_id, 
                a.lowest_taxon_id AS taxon_id,
                a.asm_name AS assembly_name,
                s.common_name, 
                a.refseq_accession AS assembly_refseq_accession,
                a.release_date AS assembly_date,
                s.scientific_name AS species_name,
                a.assembly_id, 
                mb.bioproject_name AS assembly_group,
                CONCAT(a.gca_chain, '.', a.gca_version) AS gca
            FROM assembly a
            JOIN bioproject b ON a.assembly_id = b.assembly_id
            JOIN species s ON a.lowest_taxon_id = s.lowest_taxon_id
            LEFT JOIN main_bioproject mb ON b.bioproject_id = mb.bioproject_id
            WHERE CONCAT(a.gca_chain, '.', a.gca_version) IN ({placeholders})
        """

        registry_info = mysql_fetch_data(
            registry_query,
            host=server_info["registry"]["db_host"],
            user=server_info["registry"]["db_user"],
            port=server_info["registry"]["db_port"],
            database="gb_assembly_metadata",
            password="",
            params=assembly_accessions,
        )

        if registry_info:
            return registry_info[0]
        return None

    except pymysql.Error as err:
        logger.error(f"MySQL error: {err}")
        return {}


def add_generated_data(server_info, assembly_accession, init_file, clade_settings_path, production_gca):
    registry_info = get_metadata_from_registry(server_info, assembly_accession, init_file)
    if registry_info is None:
        raise ValueError("registry_info is None")

    clade, genus_id, clade_metadata = assign_clade(
        server_info,
        registry_info,
        clade_settings_path
    )

    registry_info["clade"] = clade
    registry_info["genus_taxon_id"] = genus_id



    if clade_metadata:
        registry_info.update(clade_metadata)

    info_dict = registry_info

    # Create variables for pipeline
    info_dict["strain_type"] = "strain"

    if "alternate_haplotype" in registry_info.get("assembly_name", ""):
        info_dict["common_name"] = "alternate haplotype"
        info_dict["species_strain"] = "alternate haplotype"

    if not registry_info.get("common_name"):
        info_dict["common_name"] = "NA"

    info_dict["species_url"] = (
        f"{registry_info['species_name']}_{assembly_accession}"
    )
    info_dict["species_display_name"] = (
        f"{registry_info['species_name']} ({registry_info['common_name']}) - {assembly_accession}"
    )
    info_dict["species_strain"] = "reference"

    raw_species = (
        registry_info.get("species_name", "").strip().lower().replace(" ", "_")
    )
    species_name = raw_species.rstrip("_")  # Remove trailing underscore

    # Extract binomial name
    parts = species_name.split("_")
    if len(parts) >= 2:
        p1, p2 = parts[:2]
        binomial_species_name = f"{p1}_{p2}"
        max_len = 15
        production_name = f"{p1[:max_len]}_{p2[:max_len]}"
    else:
        binomial_species_name = ""
        production_name = ""

    production_name += f"_{production_gca}"

    # Update dictionary
    info_dict["species_name"] = species_name
    info_dict["binomial_species_name"] = binomial_species_name
    info_dict["production_name"] = production_name
    info_dict["species_strain_group"] = production_name

    return info_dict



def custom_loading(init_file):

    # Initialize variables
    custom_dict = {}

    if init_file.exists():
        try:
            with init_file.open("r") as f:
                logger.info("Using custom loading .ini file.")

                for line in f:
                    line = line.strip()

                    if "=" in line:
                        key, value = map(str.strip, line.split("=", 1))
                        logger.info(f"Found key/value pair: {key} => {value}")
                        custom_dict[key] = value
                    elif line == "":
                        continue
                    else:
                        logger.error(f"Line format not recognised. Skipping line:\n{line}")

        except OSError as e:
            raise Exception(f"Could not open or read {init_file}") from e

    return custom_dict


def create_dir(path, mode=None):
    try:
        os.makedirs(path, exist_ok=True)
        if mode is not None:
            os.chmod(path, mode)
    except Exception as e:
        raise RuntimeError(f"Failed to create dir: {path}") from e
