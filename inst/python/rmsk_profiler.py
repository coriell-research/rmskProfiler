"""
Gennaro Calendo
2024-08-24

Functions used in the rmskProfiler R package for generating Salmon indeces of 
unique RepeatMasker elements
"""
import sys
import hashlib
import json
import gzip

from pathlib import Path
from collections import namedtuple
from collections import defaultdict

import pybedtools


def rmsk2bed(rmsk, exclude, min_len):
  """Convert the RepeatMasker outfile to a BED file. There are already other 
  utilities and BED files available from UCSC that exist for mouse and human 
  however, this function allows for pre-filtering of only certain elements to be
  included in the resulting BED file which simplifies extraction to a fasta 
  later
  """
  rmsk = Path(rmsk)
  outfile = rmsk.with_suffix(".bed")
  
  keep_chroms = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                 "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                 "chrX", "chrY"}
  
  f = gzip.open(rmsk, "rt")
  out = open(outfile, "w")
  
  # Skip the header lines
  for _ in range(3):
      next(f)
      
  # Processing ------------------------------------------------------------------
  line_counter = 1
  processed_lines = 0
  bad = 0
  too_short = 0
  excluded = 0
  non_canonical = 0
  for line in f:
      line_counter += 1
      l = line.strip().split()

      # Skip bad lines
      if len(l) != 15:
          bad += 1
          print(f"Bad line found at line number: {line_counter}")
          print("The offending line is:")
          print(line)
          print("Skipping...")
          continue

      score = l[0]
      chrom = l[4]
      start = int(l[5]) - 1  # Convert to 0-based BED style starts
      end = l[6]
      strand = "-" if l[8] == "C" else "+"
      rep_elem = l[9]
      class_fam = l[10].split("/")
      rep_class = class_fam[0]
      rep_id = l[14]

      # Skip non-canconical chromosomes
      if chrom not in keep_chroms:
          non_canonical += 1
          continue

      # Fill repeat family name with class name if absent
      if len(class_fam) == 1:
          rep_fam = rep_class
      else:
          rep_fam = class_fam[1]

      # Skip lines that have element in exclusion list
      if (rep_elem in exclude) or (rep_class in exclude) or (rep_fam in exclude):
          excluded += 1
          continue

      # Skip lines where the sequence isnt long enough to hash
      if abs(int(start) - int(end)) < min_len:
          too_short += 1
          continue

      # Create a BED formatted record with a score column
      name = ".".join([rep_id, rep_class, rep_fam, rep_elem])
      out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
      processed_lines += 1

  f.close()
  out.close()

  # Summary ---------------------------------------------------------------------
  print(f"\nMalformed RepeatMasker input lines (skipped): {bad}")
  print(f"Excluded records: {excluded}")
  print(f"Records shorter than {min_len} bp: {too_short}")
  print(f"Records in non-canonical chromosomes (skipped): {non_canonical}")
  print(f"Records written to file: {processed_lines}")


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
    """Extracts unique RepeatMasker sequences by first creating a dictionary """    
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

