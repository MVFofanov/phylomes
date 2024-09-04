import logging
import os


def setup_logging(output_dir: str, cluster_name: str, logging_level=logging.INFO) -> None:
    """Set up logging configuration to save logs in the specified output directory, including the cluster name."""
    log_file_path = os.path.join(output_dir, f'{cluster_name}_log_tree_analysis.log')

    # Get the root logger
    logger = logging.getLogger()

    # Remove existing handlers to avoid duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()

    # Define the logging format including the cluster name
    # log_format = f'%(asctime)s - %(levelname)s - {cluster_name} - %(message)s'
    log_format = f'%(asctime)s\t%(levelname)s\t{cluster_name}\t%(message)s'

    # Set up logging to file with the updated format
    logging.basicConfig(
        filename=log_file_path,
        level=logging_level,
        format=log_format,
        filemode='w'  # Overwrite log file for each cluster run
    )

    # Optional: Set up console logging to print logs to the console (only once)
    console = logging.StreamHandler()
    console.setLevel(logging_level)
    formatter = logging.Formatter(log_format)
    console.setFormatter(formatter)
    logger.addHandler(console)
