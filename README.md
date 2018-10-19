# vfdbQuery

### Requirements
- Python 3.6
- [ncbi-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (makeblastdb and blastn must be in your $PATH)

### Installation
```
pip install vfdbQuery
```

### Usage
```
Usage: vfdbQuery [OPTIONS]

  vfdbQuery is a simple script for querying an input genome assembly against
  the Virulence Factor Database (VFDB).

Options:
  -i, --infile PATH     FASTA file that you want to search against VFDB
                        [required]
  -db, --database PATH  Path to Virulence Factor Database (VFDB)  [required]
  --help                Show this message and exit.

```