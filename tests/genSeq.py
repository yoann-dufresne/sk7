import random


def generate_random_sequence(filename, size=10000):
    file = open(filename, "w")
    nucleotide = ['A', 'C', 'G', 'T']
    for i in range(size):
        file.write(random.choice(nucleotide))
        if i % 60 == 0 and i > 0:
            file.write('\n')
    file.close()

