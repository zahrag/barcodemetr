
import numpy as np
import random
import os



class BarcodeSampler(object):
    def __init__(self, dna_barcodes, chunk_size):

        self.chunk_size = chunk_size

        # Get the maximum length of the sequences
        self.max_length = max(len(seq) for seq in dna_barcodes)

        # Split sequences into chunks
        self.chunked_sequences_total = [dna_barcodes[i:i + chunk_size] for i in range(0, len(dna_barcodes), chunk_size)]

        # Save selecetd chunks indices
        self.chunk_ids = [i for i in range(len(self.chunked_sequences_total))]

    def random_sampling(self, dna_barcodes):
        """
         Randomly select a chunk of sequence
        :param dna_barcodes:
        :return:
        """
        random_chunk = random.choice(self.chunk_ids)
        self.chunk_ids.remove(random_chunk)
        chunked_sequences = self.chunked_sequences_total[random_chunk]
        return chunked_sequences

