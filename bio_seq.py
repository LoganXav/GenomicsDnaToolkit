from bio_structs import CODON_LENGTH, CODON_TABLE, NUCLEOTIDE_BASE
from collections import Counter
import random

class bio_seq:

    """DNA sequence class. Handles DNA sequences and validates them."""
    def __init__(self, seq="ATCG", type="DNA", label="No label"):
        """Initialize the DNA sequence with a default sequence, type, and label."""
        self.seq = seq.upper()  # Convert sequence to uppercase
        self.label = label  # Label for the sequence
        self.type = type  # Type of sequence (DNA, RNA, etc.)
        self.is_valid = self._validate()  # Validate the sequence
        assert self.is_valid, f"Provided data does not seem to be correct {self.type}"

    def _validate(self):
        """
        Private method to check if the sequence is a valid DNA string.
        It verifies that the sequence only contains valid nucleotide bases.
        
        Returns:
            bool: True if the sequence is valid, False otherwise.
        """
        valid_bases = NUCLEOTIDE_BASE[self.type]
        return set(valid_bases).issuperset(self.seq)
    
    def get_sequence_info(self):
        """
        Returns a string containing the full sequence information.
        Output format: [Label]: <label>\n[Sequence]: <sequence>\n[Biotype]: <biotype>
        """
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.type}"
    
    def get_bio_type(self):
        """"Returns sequence type"""
        return self.type
    
    def generate_random_seq(self, length=30, type="DNA"):
        """
        Generates a random sequence of specified length and type.
        Default length is 30 and type is DNA.
        The sequence is generated using random selection of appropriate nucleotide bases.
        The class is reinitialized with the new sequence.
        """
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[type]) for char in range(length)])
        # reinitialize the class because we now have a new dna string
        self.__init__(seq, type, "Randomly generated sequence")

    def nucleotide_frequency(self):
        """
        Returns a dictionary containing the frequency count of each nucleotide in the sequence.
        Utilizes the Counter class from the collections module.
        """
        return dict(Counter(self.seq))
    
    def transcription(self):
        """DNA -> RNA Transcription. Returns a new sequence replacing Thymine nucleotides with Uracil"""
        return self.seq.replace("T", "U") if self.type == "DNA" else "Not a DNA sequence"
    
    # def complement(self):
    #     return ''.join(DNA_COMPLEMENT[char] for char in self.seq) 
    
    def reverse_complement(self):
        """Taking the complement of each nucleotide and reversing the resulting string"""
        if self.type == "DNA":
            mapping = str.maketrans('ACGT', 'TGCA')
        else:
            mapping = str.maketrans('ACGU', 'UGCA')
        return self.seq.translate(mapping)[::-1]
    
    def gc_content(self):
        """Returns the percentage of occurence of the Cytosine and Guanine nucleotide bases."""
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)
    
    def gc_content_sub_seq(self, k=20):
        """Returns the percentage of occurence of the Cytosine and Guanine nucleotide bases in sub sequences of a given k-mer length."""
        res = []
        # step size is explicitely defined as k. By default, range uses a step size of 1.
        # meaning it jumps k bases forward in each iteration
        for i in range(0, len(self.seq) - k + 1, k):
            sub_seq = self.seq[i:i+k]
            res.append(f"{ round((sub_seq.count('C') + sub_seq.count('G')) / len(sub_seq) * 100)}%")
        return res
    
    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an Amino acid sequence"""
        return [CODON_TABLE[self.type][self.seq[pos: pos + CODON_LENGTH]] for pos in range(init_pos, (len(self.seq) - CODON_LENGTH) + 1, CODON_LENGTH) ]
        
    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        # List of all the codon encodings that match the given amino acid --> ['CTG', 'CTC']
        tmp_list = [] 
        for i in range(0, (len(self.seq) - CODON_LENGTH) + 1, CODON_LENGTH):
            if CODON_TABLE[self.type][self.seq[i:i+CODON_LENGTH]] == aminoacid:
                tmp_list.append(self.seq[i:i+CODON_LENGTH])

        # Dictionary of each encoding and the number of occurences --> {'CTG': 2} 
        freq_dict = dict(Counter(tmp_list))
        total_wight = sum(freq_dict.values())
        for self.seq in freq_dict:
            # Frequency of each matched codon encoding --> {'CTG': 1}
            freq_dict[self.seq] = round(freq_dict[self.seq] / total_wight, 2)
        return freq_dict
    
    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including the reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames
    
    def __proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an amino acid sequence and return a list of possible proteins"""
        currentProtein = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # STOP accumulating amino acids
                if currentProtein: # if the current protein list isn't empty
                    for p in currentProtein:
                        proteins.append(p) # add all the proteins so far and empty the list
                    currentProtein = []
            else:
                # START accumulating amino acids
                if aa == "M":
                    # everytime it finds an M, it creates a new entry in the list and accumulates each entry until it finds a stop codon that resets it
                    currentProtein.append("")
                for i in range(len(currentProtein)):
                    currentProtein[i] += aa
        return proteins
    
    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""

        #  this step check if we want to generate reading frames from a subsequence or the entire sequence
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos:endReadPos], self.type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = [] # storing here so we can sort the list
        for rf in rfs:
            proteins = self.__proteins_from_rf(rf)
            for p in proteins:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)
        return res
