def readTextfile(filePath):
    with open(filePath, "r") as f:
        return "".join([l.strip() for l in f.readlines()])

def writeTextfile(filePath, seq, mode="w"):
    if seq is None:
        raise ValueError("The 'seq' parameter cannot be None.")
    with open(filePath, mode) as f:
        f.write(seq + "\n")

def read_FASTA(filePath):
    with open(filePath, "r") as f:
        FASTAFile = [l.strip() for l in f.readlines()]
    FASTADict = {}
    FASTALabel = ""

    # Populates the FASTADict with the sequence and its ID
    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else: 
            # PROBLEM - this assumes the FASTA file is well formatted because it just appends every sequence to the preceeding dict label until it meets another header.
            FASTADict[FASTALabel] += line
    return FASTADict
