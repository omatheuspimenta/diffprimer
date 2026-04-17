import sys
import os
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path("src").resolve()))
from diffprimer.report import generate_html_report

os.makedirs("mock_results", exist_ok=True)

# Create a rich dummy dataset
mock_data = {
    "Header": [
        "Sequence_A_region_1_100_500",
        "Sequence_A_region_2_1000_1500",
        "Sequence_B_region_1_200_600",
        "Sequence_C_region_1_50_400"
    ],
    "Forw_Seq": ["ACGTACGTACGT", "TGCATGCATGCA", "ATGCATGCATGC", "GCATGCATGCAT"],
    "Rev_Seq": ["TACGTACGTACG", "CATGCATGCATG", "ATGCATGCATGC", "TGCATGCATGCA"],
    "Amplicon_Seq": ["ACGTACGTACGTNNNNNNNNNNNNNCGTACGTACGTA", "TGCATGCATGCANNNNNNNNNNNNNNCATGCATGCATG", "ATGCATGCATGCNNNNNNNNNNNNNNGCATGCATGCAT", "GCATGCATGCATNNNNNNNNNNNNNTGCATGCATGCA"],
    "Amplicon_Size": [37, 38, 38, 37],
    "Region_Length": [400, 500, 400, 350],
    "Annotation": [
        "ID=gene001;Name=kin1;product=Kinase",
        "ID=gene002;product=Polymerase",
        "no annotation",
        "gene=hypothetical_003"
    ],
    "Sequence_Region": ["NNNACGTACGTACGTNNNNNNNNNNNNNCGTACGTACGTANNN", "NNNTGCATGCATGCANNNNNNNNNNNNNNCATGCATGCATGNNN", "NNNATGCATGCATGCNNNNNNNNNNNNNNGCATGCATGCATNNN", "NNNGCATGCATGCATNNNNNNNNNNNNNTGCATGCATGCANNN"],
    "Specificity_Tag": [
        "Specific_LowGlobalSim",
        "Specific_LowGlobalSim",
        "Specific_PositionalMismatches",
        "NonSpecific_Amplification"
    ],
    "Most_Similar_Target": ["-", "-", "OffTarget_Seq_B", "OffTarget_Seq_C"],
    "Max_Similarity_pct": ["-", "-", "78.5", "99.0"],
    "Left_TM": [60.5, 61.2, 59.8, 62.0],
    "Right_TM": [60.1, 60.8, 60.4, 61.5],
    "Left_GC": [50.0, 55.0, 45.0, 60.0],
    "Right_GC": [52.0, 50.0, 48.0, 58.0],
    "Left_HAIRPIN_TH": [30.0, 40.0, 25.0, 48.0],
    "Right_HAIRPIN_TH": [32.0, 38.0, 28.0, 45.0],
    "Left_SELF_ANY_TH": [20.0, 25.0, 18.0, 35.0],
    "Right_SELF_ANY_TH": [22.0, 20.0, 15.0, 30.0],
}

df = pd.DataFrame(mock_data)
df.to_csv("mock_results/primer_design_results.csv", sep=";", index=False)

# Contig info
contig_info = [
    {
        "header": "Sequence_A",
        "length": 5000,
        "regions": [
            {"start": 100, "end": 500},
            {"start": 1000, "end": 1500}
        ]
    },
    {
        "header": "Sequence_B",
        "length": 8000,
        "regions": [
            {"start": 200, "end": 600}
        ]
    },
    {
        "header": "Sequence_C",
        "length": 4000,
        "regions": [
            {"start": 50, "end": 400}
        ]
    },
    {
        "header": "Sequence_D_NoPrimers",
        "length": 3000,
        "regions": [
            {"start": 10, "end": 200}
        ]
    }
]

run_params = {
    "Reference File": "/long/path/to/my/reference.fasta",
    "Sequences Path": "/long/path/to/my/target_sequences/",
    "Annotation File": "/long/path/to/my/annotations.gff3",
    "K-mer Size": "21",
    "Min Region Length": "200"
}

generate_html_report(
    results_dir="mock_results",
    run_params=run_params,
    primer3_args={"PRIMER_OPT_SIZE": "20", "PRIMER_MIN_TM": "57.0"},
    contig_info=contig_info,
    version="1.0.0-mock"
)
print("Mock generated.")
