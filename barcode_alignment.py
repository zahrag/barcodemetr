from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import os
from pathlib import Path

# Get the path of the currently running script
current_directory = Path(__file__).parent

def make_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def align_with_mafft(sequences):
    """
    This function performs MAFFT alignment of a chunk of sequences without saving to disk.
    :param sequences: Chunk of sequences.
    :param taxa: Taxonomic rank.
    :param max_length: Maximum length of sequences.
    :return: Aligned sequences.
    """

    # Get the maximum length of the sequences
    max_length = max(len(seq) for seq in sequences)

    print(f"Aligning chunk of sequences with {len(sequences)} sequences.")

    # Format sequences into FASTA format in-memory
    fasta_io = StringIO()
    for j, seq in enumerate(sequences):
        # Optionally trim or pad sequences here
        if len(seq) > max_length:
            seq = seq[:max_length]  # Trim to max length
        else:
            seq = seq + '-' * (max_length - len(seq))  # Pad to max length
        fasta_io.write(f">sequence_{j}\n{seq}\n")
    
    # Reset the buffer's position to the beginning
    fasta_io.seek(0)

    # Run MAFFT alignment
    mafft_cline = MafftCommandline(input="")
    stdout, stderr = mafft_cline(stdin=fasta_io.getvalue())
    
    # Read the alignment from stdout directly
    aligned_io = StringIO(stdout)
    aligned_sequences = AlignIO.read(aligned_io, "fasta")
    aligned_sequences = [str(record.seq) for record in aligned_sequences]

    return aligned_sequences


def mafft_align_and_save(sequences, taxa, remove_fasta=False):
    """
    This function performs MAFFT alignment of a chunk of sequences.
    :param sequences: A chunk of sequences.
    :param taxa: Taxonomic rank.
    :return: Aligned sequences.
    """
    # Get the maximum length of the sequences
    max_length = max(len(seq) for seq in sequences)

    # print(f"Aligning chunk of sequences with {len(sequences)} sequences.")
    fasta_dir = f'{current_directory}/fasta_files/{taxa}'
    make_directory(fasta_dir)
    fasta_filename = f"{fasta_dir}/{taxa}_sequences_chunk.fasta"
    aligned_filename = f"{fasta_dir}/{taxa}_aligned_chunk.fasta"

    # Write the chunk of sequences to a FASTA file
    with open(fasta_filename, "w") as file:
        for j, seq in enumerate(sequences):
            # Optionally trim or pad sequences here
            if len(seq) > max_length:
                seq = seq[:max_length]  # Trim to max length
            else:
                seq = seq + '-' * (max_length - len(seq))  # Pad to max length
            file.write(f">sequence_{j}\n{seq}\n")

    # Run MAFFT alignment
    mafft_cline = MafftCommandline(input=fasta_filename)
    stdout, stderr = mafft_cline()

    # Save the aligned output
    with open(aligned_filename, "w") as aligned_file:
        aligned_file.write(stdout)

    # Read aligned sequences for further processing
    aligned_sequences = AlignIO.read(aligned_filename, "fasta")
    aligned_sequences = [str(record.seq) for record in aligned_sequences]

    if remove_fasta:
        os.remove(fasta_filename)
        os.remove(aligned_filename)

    return aligned_sequences


def perform_mafft_alignment(sequences, taxa, save_fasta=True):

    if save_fasta:
        return mafft_align_and_save(sequences, taxa,  remove_fasta=False)

    align_with_mafft(sequences)





