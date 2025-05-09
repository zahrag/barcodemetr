from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import os


class BarcodeAlignment:
    def __init__(self, max_seq_length=None, save_path=None, remove_fasta=False):
        self.save_path = save_path
        self.max_seq_length = max_seq_length  # Applied when using a fixed length for all barcodes
        self.remove_fasta = remove_fasta   # Applied when lacking space to keep fasta files

    def mafft_align(self, sequences):
        """
        This function performs MAFFT alignment of a chunk of sequences without saving fasta files to disk.
        :param sequences: Chunk of sequences.
        :return: Aligned sequences.
        """

        # Get the maximum length of the sequences
        max_length = (
        self.max_seq_length if self.max_seq_length is not None 
        else max(len(seq) for seq in sequences)
        )
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


    def mafft_align_and_save(self, sequences, taxa):
        """
        This function performs MAFFT alignment of a chunk of DNA barcode sequences saving fasta files on disk.
        :param sequences: A chunk of sequences.
        :param taxa: Taxonomic rank.
        :return: Aligned sequences.
        """
        # Get the maximum length of the sequences
        max_length = (
        self.max_seq_length if self.max_seq_length is not None 
        else max(len(seq) for seq in sequences)
        )

        # print(f"Aligning chunk of sequences with {len(sequences)} sequences.")
        fasta_dir = f'{self.save_path}/fasta_files/{taxa}'
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)

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

        if self.remove_fasta:
            os.remove(fasta_filename)
            os.remove(aligned_filename)

        return aligned_sequences


    def perform_mafft_alignment(self, sequences, taxa):

        if self.save_path is None:
            return self.mafft_align(sequences)

        return self.mafft_align_and_save(sequences, taxa)






