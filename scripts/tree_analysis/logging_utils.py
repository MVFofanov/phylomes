import logging
import matplotlib
import os


def setup_logging(output_dir: str, logging_level=logging.INFO) -> None:
    """Set up logging configuration to save logs in the specified output directory."""
    log_file_path = os.path.join(output_dir, 'log_tree_analysis.log')
    logging.basicConfig(
        filename=log_file_path,
        level=logging_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filemode='w'
    )

# Set logging level for matplotlib to WARNING or higher
logging.getLogger('matplotlib').setLevel(logging.WARNING)
