import unittest
import random
import subprocess
import genSeq
import time


class MyTestCase(unittest.TestCase):

    def test_kmer_value_and_mini(self):
        def find_value(kmer):
            res = 0
            encoding = {"A": 0, "C": 1, "G": 2, "T": 3}
            for c in kmer:
                res = (res << 2) + encoding[c]
            return res

        def find_mini(kmer, size):
            minimiser = "T" * size
            for i in range(len(kmer) - size + 1):
                if kmer[i:i + size] < minimiser:
                    minimiser = kmer[i:i+size]
            return minimiser

        def parse_out(filename):
            file = open(filename, "r")
            file.readline() # skip the first line
            for line in file:
                listing = line.split()
                kmer = listing[2]
                value = int(listing[7])
                minimiser = listing[-6]
                mini_val = int(listing[-1])
                self.assertEqual(find_value(kmer), value)
                self.assertEqual(find_mini(kmer, len(minimiser)), minimiser)
                self.assertEqual(find_value(minimiser), mini_val)
            file.close()

        genSeq.generate_random_sequence('tests/seq.txt', 1000000)
        print("Sequence generated")

        k = random.randint(15, 25)
        m = random.randint(5, 10)
        inputs = ['tests/seq.txt', str(k), str(m)]
        cmd = ["./bin/testKmerMini"] + inputs
        out = open("tests/outfile.txt", "w")
        start_time = time.time()
        subprocess.run(cmd, stdout=out)
        exec_time = time.time() - start_time
        print(f"Got cpp output in {exec_time}s")
        out.close()
        parse_out("tests/outfile.txt")


if __name__ == '__main__':
    unittest.main()
