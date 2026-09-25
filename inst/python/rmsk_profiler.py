"""
Gennaro Calendo
2024-08-24

Functions used in the rmskProfiler R package for generating Salmon indeces of 
unique RepeatMasker elements
"""
import hashlib
import json

from pathlib import Path
from collections import namedtuple
from collections import defaultdict

import pybedtools


def read_fasta(fa):
    """Read records from a fasta file"""
    FastaRecord = namedtuple("FastaRecord", ["header", "sequence"])
    header = ""
    sequence = ""
    with open(fa, "r") as f:
        for line in f:
            line = line.strip()
            if line[0] == ">":
                if sequence:
                    yield FastaRecord(header, sequence)
                header = line[1:]
                sequence = ""
            else:
                sequence += line
        yield FastaRecord(header, sequence)


def extract_unique_records(bed, fasta):
    """Extracts unique RepeatMasker sequences by first extracting all sequences from
    the RepeatMasker BED file into a fasta and then collecting the unique sequences
    by SHA1 hash. These unique sequences are then dumped to a json file and a fasta
    file.
    """    
    bedfile = Path(bed)
    fasta = Path(fasta)
    out_json = Path(bedfile.parent, "rmsk-duplicateInfo.json")
    out_fa = Path(bedfile.parent, "rmsk-unique.fa")
    
    print("Reading in the BED file...")
    a = pybedtools.BedTool(bedfile)
    print("Extracting RepeatMasker sequences...")
    a = a.sequence(fi=str(fasta), name=True, s=True)
    
    print("Hashing RepeatMasker sequences and checking for duplicates...")
    hash_dict = defaultdict(list)
    for record in read_fasta(a.seqfn):
        seq = record.sequence
        header = record.header
        seq_b = seq.encode()
        m = hashlib.sha1()
        m.update(seq_b)
        h = m.hexdigest()
        hash_dict[h].append(header)

    print("Dumping duplicate information to json...")
    with open(out_json, "w") as out:
        json.dump(hash_dict, out, indent="\t")

    print("Getting unique fasta records...")
    uniq_headers = [l[0] for l in hash_dict.values()]
    uniq_hashes = hash_dict.keys()
    uniq_dict = dict(zip(uniq_headers, uniq_hashes))

    print("Creating fasta file of unique records...")
    keep_headers = []
    keep_sequences = []
    for record in read_fasta(a.seqfn):
        header = record.header
        seq = record.sequence
        if header in uniq_dict.keys():
            hash_val = uniq_dict[header]
            h = ">" + hash_val
            s = seq
            keep_headers.append(h)
            keep_sequences.append(s)

    print("Writing out unique fasta records to file...")
    records = list(zip(keep_headers, keep_sequences))
    with open(out_fa, "w") as out2:
        out2.write("\n".join("{}\n{}".format(x[0], x[1]) for x in records))

    print("Done.")

