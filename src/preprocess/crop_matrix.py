"""
Crop extra elements to the right of a NBDB matrix file
"""

def crop(input_file, output_file, count=20):
    matrix = []
    with open(input_file, 'r') as file:
        for line in file:
            matrix.append(line.strip().split(" ")[:count])

    with open(output_file, 'w') as file:
        for line in matrix:
            file.write(" ".join(line) + "\n")
