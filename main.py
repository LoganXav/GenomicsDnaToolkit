from bio_seq import bio_seq
from utils import writeTextfile, readTextfile, read_FASTA

# Create an instance of the bio_seq class
my_seq = bio_seq(seq="ATGCGTAGCTAGCTAGCTAGCTAG", type="DNA", label="Test Sequence")

# Call the public methods and print their outputs
print(my_seq.get_sequence_info())
print(my_seq.get_bio_type())
my_seq.generate_random_seq(length=50, type="RNA")
print(my_seq.get_sequence_info())
print(my_seq.nucleotide_frequency())
print(my_seq.transcription())
print(my_seq.reverse_complement())
print(my_seq.gc_content())
print(my_seq.gc_content_sub_seq(k=5))
print(my_seq.translate_seq(init_pos=0))
print(my_seq.codon_usage(aminoacid="L"))
print(my_seq.gen_reading_frames())
print(my_seq.all_proteins_from_orfs(startReadPos=0, endReadPos=0, ordered=False))

my_seq.generate_random_seq(50)
writeTextfile("test.txt", my_seq.seq)
for rf in my_seq.gen_reading_frames():
    writeTextfile("test.txt", str(rf), "a")

fa = read_FASTA("sample_fasta.txt")

print(fa)
