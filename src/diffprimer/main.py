from diffprimer.kmer_extractor import process_seqs, PrimerCandidate
from diffprimer.config import load_config
from diffprimer.helpers import write_csv, get_primers, write_csv_header, annotation_dataframe, get_sequence_intersection
from diffprimer.logs import diffprimerLog, _make_progress
import os
import multiprocessing
from functools import partial
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.columns import Columns
from rich import box


import warnings

# Suppress harmless RuntimeWarning from PyO3/Python 3.14+ free-threading compatibility checks
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*global interpreter lock.*")

logger = diffprimerLog()
console = Console()

   
# Wrapper for parallel processing
def design_primers_for_contig(contig_data, global_args):
    """
    Worker function to design primers for a single contig.
    Returns a list of result dictionaries to be written to CSV.
    """
    header, sequence, regions = contig_data
    results = []
    
    n_region = 0
    for region in regions:
        # region is passed as a dict to avoid pickling issues with PyO3 objects
        start = region["start"]
        end = region["end"]
        
        region_seq = sequence[start:end]
        get_primers_result = get_primers(region_seq, global_args)
        
        if (
            get_primers_result["best_left_primer"]
            and get_primers_result["best_right_primer"]
        ):
            # Prepare data for writing
            # We don't write here to avoid locking issues
            try:
                # We need to compute the specific header for this region
                region_header = f"{header}_region_{n_region + 1}_{start}_{end}"
                results.append({
                    "result_dict": get_primers_result,
                    "header": region_header,
                    "sequence": region_seq,
                })
            except Exception as e:
                logger.error(f"Error processing region in {header}: {e}")
                
    return results

def main(
    reference_file: str,
    sequences_path: str,
    annotation_path: str,
    config_file: str,
    k: int = 21,
    cpus: int | None = None,
    min_region_length: int = 200,
    reference_max_abundance: int = 1,
    check_specificity: bool = False,
    similarity_threshold: float = 80.0,
) -> None:    
    input_table = Table(show_header=False, box=None, padding=(0, 2))
    input_table.add_column("Parameter", style="bold cyan")
    input_table.add_column("Value", style="white")

    input_table.add_row("Reference File", str(reference_file))
    input_table.add_row("Sequences Path", str(sequences_path))
    input_table.add_row("Annotation File", str(annotation_path))
    input_table.add_row("Configuration File", str(config_file))
    input_table.add_row("k-mer Size", str(k))
    input_table.add_row("Min Region Length", str(min_region_length))
    input_table.add_row("Ref Max Abundance", str(reference_max_abundance))
    input_table.add_row("Check Specificity", str(check_specificity))

    if check_specificity:
        input_table.add_row("Similarity Threshold", str(similarity_threshold))

    console.print(Panel(
        input_table, 
        title="[bold green]DiffPrimer Settings[/bold green]", 
        border_style="green",
        expand=False
    ))
    
    # --- 2. Create a 4-Column Table for Primer3 Args ---
    # We use 4 columns to create the visual effect of "Two Columns of Data"
    p3_table = Table(
        show_header=True, 
        header_style="bold blue", 
        box=box.SIMPLE_HEAD,
        expand=True
    )
    p3_table.add_column("Parameter", style="dim")
    p3_table.add_column("Value", style="yellow")
    p3_table.add_column("Parameter", style="dim")
    p3_table.add_column("Value", style="yellow")
    
    global_args = load_config(config_file)
    # Logic to split dictionary into pairs for the table
    items = list(global_args.items())
    for i in range(0, len(items), 2):
        key1, val1 = items[i]
        # Check if there is a second item in this pair
        if i + 1 < len(items):
            key2, val2 = items[i+1]
            p3_table.add_row(key1, str(val1), key2, str(val2))
        else:
            # If odd number of items, leave the right side empty
            p3_table.add_row(key1, str(val1), "", "")

    # --- 3. Print Everything ---
    

    console.print(Panel(
        p3_table, 
        title="[bold blue]Primer3 Global Arguments[/bold blue]", 
        border_style="blue",
        expand=False
    ))
    
    if cpus is None:
        detected = os.cpu_count()
        if detected is None:
            raise RuntimeError(logger.error("Could not detect number of CPUs available. Please specify --cpus manually."))
        cpus = detected - 1
    logger.info(f"Using {cpus} CPU cores for parallel processing.")
    
    input_file = os.path.abspath(sequences_path)
    results_dir = os.path.join(os.path.dirname(input_file), "results")
    
    # Load annotation dataframe once
    logger.info("Loading annotation file...")
    try:
        df_annotation = annotation_dataframe(annotation_path)
    except Exception as e:
         logger.warning(f"Could not load annotation file: {e}")
         df_annotation = None

    os.makedirs(results_dir, exist_ok=True)

    output_file = os.path.join(results_dir, "primer_design_results.csv")
    write_csv_header(output_file)
    
    regions_output_path = os.path.join(results_dir, "exclusive_regions.txt")
    fasta_output_path = os.path.join(results_dir, "exclusive_sequences.fasta")
    
    # process_seqs now handles its own progress with indicatif (Rust)
    # We don't wrap it in a Python spinner to avoid cursor conflicts
    result = process_seqs(
        reference_path=reference_file,
        sequences_dir=sequences_path,
        k=k,
        num_threads=cpus,
        min_region_length=min_region_length,
        regions_output_path=regions_output_path,
        fasta_output_path=fasta_output_path,
        max_abundance=reference_max_abundance,
    )
    
    # Prepare data for multiprocessing
    # Extract simple data structures to avoid pickling issues with PyO3 objects if any
    contig_data_list = []
    for contig in result:
        # Convert Rust Region objects to pure Python types (dicts)
        regions_data = [{"start": r.start, "end": r.end} for r in contig.regions]
        contig_data_list.append((contig.header, contig.sequence, regions_data))

    logger.info(f"Designing primers for {len(contig_data_list)} contigs using {cpus} workers...")

    # Collect all designed primers first
    designed_primers_list = []

    with multiprocessing.Pool(processes=cpus) as pool:
        # Use partial to pass constant global_args
        worker_func = partial(design_primers_for_contig, global_args=global_args)
        
        with _make_progress(total=True) as progress:
            task = progress.add_task("[green]Designing primers...", total=len(contig_data_list))
            
            # Use imap to get results as they complete
            for batch_results in pool.imap_unordered(worker_func, contig_data_list):
                for item in batch_results:
                    designed_primers_list.append(item)
                progress.advance(task)

    # NEW: Specificity Check
    spec_map = {}
    if check_specificity and designed_primers_list:
        logger.info(f"Checking specificity for {len(designed_primers_list)} primer candidates...")
        
        candidates = []
        for item in designed_primers_list:
            res_dict = item["result_dict"]
            header = item["header"]
            seq = item["sequence"]
            
            left = res_dict.get("best_left_primer", {})
            right = res_dict.get("best_right_primer", {})
            
            left_seq = left.get("SEQUENCE", "")
            right_seq = right.get("SEQUENCE", "")
            
            # Get offsets (0-based start index)
            l_offset = left.get("COORDS", [0, 0])[0]
            r_offset = right.get("COORDS", [0, 0])[0]
            
            # Parse region start/end from header if possible
            # Format: "{original}_region_{n}_{start}_{end}"
            try:
                parts = header.split('_')
                r_end = int(parts[-1])
                r_start = int(parts[-2])
            except (ValueError, IndexError):
                r_start = 0
                r_end = len(seq)

            # Only create candidate if primers exist
            if left_seq and right_seq:
                cand = PrimerCandidate(
                    header=header,
                    region_sequence=seq,
                    start=r_start,
                    end=r_end,
                    left_primer_seq=left_seq,
                    right_primer_seq=right_seq,
                    left_primer_offset=l_offset,
                    right_primer_offset=r_offset
                )
                candidates.append(cand)
        
        # Run Rust Specificity Check
        # Default similarity threshold 80.0 if not specified
        sim_threshold = similarity_threshold 
        
        try:
            from diffprimer.kmer_extractor import check_specificity
            specificity_results = check_specificity(
                candidates=candidates,
                sequences_dir=sequences_path,
                similarity_threshold=sim_threshold,
                local_mismatch_threshold=3,
                num_threads=cpus
            )
            
            # Create lookup map
            spec_map = {res.region_header: res for res in specificity_results}
            logger.info("Specificity analysis completed successfully.")
            
        except Exception as e:
            logger.error(f"Specificity check failed: {e}")
            spec_map = {}

    elif not designed_primers_list:
        logger.warning("No primers designed.")
        spec_map = {}

    # Write Results to CSV
    logger.info("Writing results to CSV...")
    with _make_progress(total=True) as progress:
        task = progress.add_task("[green]Writing results to CSV...", total=len(designed_primers_list))
        for item in designed_primers_list:
            header = item["header"]
            tag = "Not_Checked"
            if header in spec_map:
                tag = str(spec_map[header].tag)
                
            write_csv(
                result_dict=item["result_dict"],
                header=header,
                sequence=item["sequence"],
                df_annotation=df_annotation,
                output_file=output_file,
                specificity_tag=tag
            )
            progress.advance(task)

    logger.info(f"Primer design completed. Results saved to {output_file}")
