# vfdbQuery

### Requirements
- Python 3.6
- [ncbi-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (makeblastdb and blastn must be in your $PATH)

### Installation
```
pip install vfdbQuery
```

### Usage
```bash
Usage: vfdbQuery [OPTIONS]

Options:
  -i, --infile PATH     FASTA file that you want to search against VDB
                        [required]
  -db, --database PATH  Path to Virulence Factor Database (VDB)  [required]
  --help                Show this message and exit.
```