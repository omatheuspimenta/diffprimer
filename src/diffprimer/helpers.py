"""
Helper functions for primer design and annotation processing.
"""

import primer3
import pandas as pd
import gffpandas.gffpandas as gffpd
import os
from diffprimer.logs import diffprimerLog

logger = diffprimerLog()

def get_primers(sequence: str, global_args: dict) -> dict:
    """
    Get the primers for a given sequence using Primer3 with the provided global arguments.
    Args:
        sequence (str): The DNA sequence for which to design primers.
        global_args (dict): A dictionary of global Primer3 configuration parameters.
    Returns:
        dict: A dictionary containing the best left and right primers, amplicon size, \
        and amplicon sequence.
    """
    if not sequence or not isinstance(sequence, str):
        logger.error("Invalid or empty sequence provided.")
        return {
            "best_left_primer": {}, "best_right_primer": {},
            "amplicon_size": None, "amplicon_seq": "",
        }
    
    
    sequence = sequence.upper()
    seq_args = {
        "SEQUENCE_ID": "region_of_interest",
        "SEQUENCE_TEMPLATE": sequence,
    }

    try:
        design_result_dict = primer3.design_primers(
            seq_args=seq_args,
            global_args=global_args,
        )
    except Exception as e:
        logger.error(f"Primer3 design failed for sequence: {e}")
        design_result_dict = {}
        
    # Return empty results if no primers were found
    num_returned = design_result_dict.get("PRIMER_PAIR_NUM_RETURNED", 0)
    if num_returned == 0:
        return {
            "best_left_primer": {},
            "best_right_primer": {},
            "amplicon_size": None,
            "amplicon_seq": "",
        }

    # The best primers are the first one in the returned list based in primer3 penalities.
    left_coords = design_result_dict.get("PRIMER_LEFT_0")   
    right_coords = design_result_dict.get("PRIMER_RIGHT_0") 

    if not left_coords or not right_coords:
        logger.error("Primer3 returned >0 pairs but missing coordinates.")
        return {
            "best_left_primer": {},
            "best_right_primer": {},
            "amplicon_size": None,
            "amplicon_seq": "",
        }

    L_start = left_coords[0]
    R_end = right_coords[0]

    amplicon_size = design_result_dict.get("PRIMER_PAIR_0_PRODUCT_SIZE")
    amplicon_seq = sequence[L_start : R_end + 1]
    
    best_left = {
        "COORDS": left_coords,
        "SEQUENCE": design_result_dict.get("PRIMER_LEFT_0_SEQUENCE"),
        "TM": design_result_dict.get("PRIMER_LEFT_0_TM"),
        "PENALTY": design_result_dict.get("PRIMER_LEFT_0_PENALTY"),
        "GC_PERCENT": design_result_dict.get("PRIMER_LEFT_0_GC_PERCENT")
    }

    best_right = {
        "COORDS": right_coords,
        "SEQUENCE": design_result_dict.get("PRIMER_RIGHT_0_SEQUENCE"),
        "TM": design_result_dict.get("PRIMER_RIGHT_0_TM"),
        "PENALTY": design_result_dict.get("PRIMER_RIGHT_0_PENALTY"),
        "GC_PERCENT": design_result_dict.get("PRIMER_RIGHT_0_GC_PERCENT")
    }

    return {
        "best_left_primer": best_left,
        "best_right_primer": best_right,
        "amplicon_size": amplicon_size,
        "amplicon_seq": amplicon_seq,
    }


def annotation_dataframe(
    annotation_path: str | None,
) -> pd.DataFrame | None:
    """
    Process annotation information from a GFF3 file into a DataFrame and a sequence \
    interval Series.

    Args:
        annotation_path (str): The path to the GFF3 annotation file.

    Returns:
        df_annotation (pd.DataFrame): A DataFrame with processed annotation data.
    """
    if annotation_path is None:
        logger.info("No annotation file provided. Skipping annotation processing.")
        return None  # Return None if no annotation file is provided
    df_annotation = gffpd.read_gff3(annotation_path).attributes_to_columns()

    df_annotation["interval"] = df_annotation.apply(
        lambda row: pd.Interval(row.start, row.end, closed="both"), axis=1
    )

    return df_annotation


def get_sequence_intersection(df_annotation: pd.DataFrame, header_complete: str) -> str:
    """ 
    Generate the annotation string for a given sequence region based on overlaps with \
    annotated features.
    Args:
        df_annotation (pd.DataFrame): DataFrame containing annotation data with intervals.
        header_complete (str): The complete header string containing sequence ID and \
        region coordinates.
    Returns:
        str: A string representing the overlapping annotations or 'no annotation' if \
        none are found.
    """
    row_split = header_complete.split("_")
    header = row_split[0]
    start = row_split[-2]
    end = row_split[-1]

    q_interval = pd.Interval(int(start), int(end), closed="both")
    overlaps = df_annotation[
        (df_annotation["seq_id"] == header)
        & (df_annotation["interval"].apply(lambda iv: iv.overlaps(q_interval)))
    ]
    if overlaps.empty:
        return "no annotation"
    else:
        return "|".join(list(overlaps["attributes"])).replace(";", ",")

def write_csv_header(output_file: str) -> None:
    """
    Write the header line to the output CSV file.
    Args:
        output_file (str): Path to the output CSV file.
    Returns:
        None
    """
    try:
        with open(output_file, mode="a", encoding="utf-8") as csvfile:
            # Write header line if the file is new
            file_exists = os.path.exists(output_file)
            if not file_exists or os.path.getsize(output_file) == 0:
                csvfile.write(
                    "Header;Forw_Seq;Rev_Seq;Amplicon_Seq;Amplicon_Size;"
                    "Left_Penalty;Right_Penalty;"
                    "Left_TM;Right_TM;"
                    "Left_GC;Right_GC;"
                    "Left_SELF_ANY_TH;Right_SELF_ANY_TH;"
                    "Left_SELF_END_TH;Right_SELF_END_TH;"
                    "Left_HAIRPIN_TH;Right_HAIRPIN_TH;"
                    "Left_END_STABILITY;Right_END_STABILITY;"
                    "Region_Length;"
                    "Annotation;"
                    "Sequence_Region;"
                    "Specificity_Tag\n"
                )
    except Exception as e:
        raise RuntimeError(logger.error(f"Error writing CSV header: {e}"))

def write_csv(
    result_dict: dict,
    header: str,
    sequence: str,
    df_annotation: pd.DataFrame | None, # Accept loaded DF instead of path
    output_file: str,
    specificity_tag: str = "Not_Checked",
) -> None:
    """
    Write primer design results to a CSV file.
    Args:
        result_dict (dict): Dictionary containing primer design results.
        header (str): Header string for the sequence region.
        sequence (str): The DNA sequence of the region.
        df_annotation (pd.DataFrame): Dataframe for overlap checking.
        output_file (str): Path to the output CSV file.
        specificity_tag (str): Specificity classification tag.
    Returns:
        None    
    """
    left = result_dict.get("best_left_primer", {})
    right = result_dict.get("best_right_primer", {})
    amplicon_seq = result_dict.get("amplicon_seq", "")
    
    if df_annotation is not None:
        annotation = get_sequence_intersection(df_annotation, header)
    else:
        annotation = "no annotation (file not loaded)"

    with open(output_file, mode="a") as csvfile:
        csvfile.write(
            f"{header};"
            f"{left.get('SEQUENCE', 'NA')};"
            f"{right.get('SEQUENCE', 'NA')};"
            f"{amplicon_seq};"
            f"{result_dict.get('amplicon_size', 'NA')};"
            f"{left.get('PENALTY', 'NA')};{right.get('PENALTY', 'NA')};"
            f"{left.get('TM', 'NA')};{right.get('TM', 'NA')};"
            f"{left.get('GC_PERCENT', 'NA')};{right.get('GC_PERCENT', 'NA')};"
            f"{left.get('SELF_ANY_TH', 'NA')};{right.get('SELF_ANY_TH', 'NA')};"
            f"{left.get('SELF_END_TH', 'NA')};{right.get('SELF_END_TH', 'NA')};"
            f"{left.get('HAIRPIN_TH', 'NA')};{right.get('HAIRPIN_TH', 'NA')};"
            f"{left.get('END_STABILITY', 'NA')};{right.get('END_STABILITY', 'NA')};"
            f"{len(sequence)};"
            f"{annotation};"
            f"{sequence};"
            f"{specificity_tag}\n"
        )
