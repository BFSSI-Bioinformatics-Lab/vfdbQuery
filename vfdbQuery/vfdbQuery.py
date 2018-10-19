import os
import click
import logging
import subprocess
import pandas as pd
from pathlib import Path


def convert_to_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    value = Path(value)
    return value


@click.command()
@click.option("-i", "--infile", type=click.Path(exists=True), required=True,
              help='FASTA file that you want to search against VDB', callback=convert_to_path)
@click.option("-db", "--database", type=click.Path(exists=True), required=True,
              help='Path to Virulence Factor Database (VDB)', callback=convert_to_path)
def cli(infile, database):
    logging.basicConfig(
        format='\033[92m \033[1m %(asctime)s \033[0m %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info("Started VDB Query")

    blast_file = blast(infile, database)
    filtered_file = filter_blast(blast_file)
    generate_plot(filtered_file)
    logging.info("Process Complete")


def check_db_exists(database: Path):
    # Format DB
    folder_files = list(database.parent.glob("*"))
    nhr, nin, nsq = False, False, False
    for f in folder_files:
        if f.suffix == '.nhr':
            nhr = True
        elif f.suffix == '.nin':
            nin = True
        elif f.suffix == '.nsq':
            nsq = True
    if nhr and nin and nsq:
        database_exists = True
    else:
        database_exists = False
    return database_exists


def blast(infile: Path, database: Path):
    # Create db if it doesn't already exist
    db_exists = check_db_exists(database)
    if not db_exists:
        logging.info(f"Creating BLAST database with {database.name}")
        formatdb_cmd = f"makeblastdb -in {database} -dbtype nucl"
        p = subprocess.Popen(formatdb_cmd, shell=True)
        p.wait()
    else:
        logging.info("BLAST database detected")

    # Need to perform first BLAST phase, if something is found go onto the next
    logging.info(f"Running blastn on {infile}")
    out_blast = infile.with_suffix(".VFDB_BLASTn")
    blast_cmd = f"blastn -db {database} -query {infile} -outfmt '6 qseqid stitle slen " \
                f"length qstart qend sstrand pident score' > {out_blast}"

    p = subprocess.Popen(blast_cmd, shell=True)
    p.wait()
    # Check the tab-delineated VFDB.blastn to see if either of the two essential markers are present
    activator_check = False
    with open(str(out_blast), 'r') as f:
        blastn_file = f.readlines()
    for line in blastn_file:
        target = line.split('\t')[1]
        if target == "plcR Transcriptional activator" or target == "papR Signal peptide":
            activator_check = True
    if activator_check:
        active_out_name = out_blast.with_suffix(".VFDB_Active")
        os.rename(str(out_blast), str(active_out_name))
        return active_out_name
    else:
        inactive_out_name = out_blast.with_suffix(".VFDB_NOT-Active")
        os.rename(str(out_blast), str(inactive_out_name))
        logging.info('No virulence activators found. Quitting.')
        logging.info(f'See potential virulence factors in {inactive_out_name}')
        quit()


def filter_blast(infile: Path):
    with open(str(infile), 'r') as f:
        blastn_file = f.readlines()

    filtered_out_name = infile.with_suffix(".VFDB_Active_Filtered")
    with open(str(filtered_out_name), 'w') as out:
        # Will now remove any hits that are less than 70% over 70% of the VFDB gene length
        out.write("qseqid\tstitle\tslen\tlength\tqstart\tqend\tsstrand\tpident\tscore\n")
        for line in blastn_file:
            slen = line.split('\t')[2]
            length = line.split('\t')[3]
            pident = line.split('\t')[7]
            len_percent = float(length) / float(slen)
            if len_percent >= 0.7 and float(pident) >= 70.00:
                out.write(line)
    return filtered_out_name


def generate_plot(filtered_filepath: Path):
    sample = filtered_filepath.name
    sample = sample.split('.')[0]
    df = pd.read_csv(filtered_filepath, sep='\t')

    target_dict = {
        'plcR Transcriptional activator': 0,
        'papR Signal peptide': 0,
        'HblL2 BC3104 Hemolysin BL lytic component L2': 0,
        'HblL1 BC3103 Hemolysin BL lytic component L1': 0,
        'HblB BC3102 Hemolysin BL binding component precursor': 0,
        'NheA BC1809 Non-hemolytic enterotoxin lytic component L2': 0,
        'NheB BC1810 Non-hemolytic enterotoxin lytic component L1': 0,
        'NheC BC1811 Enterotoxin C': 0,
        'CytK BC1110 Cytotoxin K': 0,
        'HlyI BC5101 Perfringolysin O precursor': 0,
        "HblB' BC3101 Hemolysin BL binding component precursor": 0,
        'HlyII BC3523 Hemolysin II': 0,
        'EntFM BC1953 Enterotoxin': 0,
        'EntA BC5239 Enterotoxin - cell wall binding protein': 0,
        'EntB BC2952 Enterotoxin - cell wall binding protein': 0,
        'EntC BC0813 Enterotoxin - cell wall binding protein': 0
    }

    for target in df['stitle']:
        if target in target_dict:
            target_dict[target] += 1

    df_processed = pd.DataFrame(list(target_dict.items()), columns=['', 'Count'])

    # Picked 16 colours from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
    my_colors = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
        '#008080', '#800000', '#aaffc3', '#808000', '#000075', '#808080'
    ]

    plt = df_processed.plot(x='', y='Count', kind='barh',
                            title=sample, xticks=range(max(df_processed['Count'] + 1)),
                            figsize=(4, 6), legend=False, color=my_colors, edgecolor='black')

    fig = plt.get_figure()

    pdf_path = filtered_filepath.parent / f'{sample}_TargetCheck.pdf'
    fig.savefig(pdf_path, bbox_inches='tight')
    png_path = filtered_filepath.parent / f'{sample}_TargetCheck.png'
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    csv_path = filtered_filepath.parent / f'{sample}_TargetCheck.tsv'
    df_processed.to_csv(csv_path, sep='\t', index=None)

    logging.info(f"Created plot at {pdf_path}")
    logging.info(f"Created plot at {png_path}")
    logging.info(f"Created CSV virulence activator count data at {csv_path}")


if __name__ == "__main__":
    cli()
