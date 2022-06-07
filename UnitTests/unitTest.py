import os
import unittest
import random
import subprocess
import genSeq
import time


class MyTestCase(unittest.TestCase):

    def test_kmer_value_and_mini(self):

        encoding = {"A": 0, "C": 1, "G": 3, "T": 2}
        def find_value(kmer):
            res = 0
            for c in kmer:
                res = (res << 2) + encoding[c]
            return res

        def lesser_than(str1:str, str2:str):
            for c1, c2 in zip(str1, str2):
                if encoding[c1] < encoding[c2]:
                    return True
                elif encoding[c1] > encoding[c2]:
                    return False
                else:
                    continue
            return True

        def find_mini(kmer, size):
            minimiser = "G" * size
            for i in range(len(kmer) - size + 1):
                if lesser_than(kmer[i:i + size],  minimiser):
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

        print("======= TESTING Minimiser =======")
        genSeq.generate_random_sequence('UnitTests/seq.txt', 1000000)
        print("Sequence generated")

        k = random.randint(15, 25)
        m = random.randint(5, 10)
        inputs = ['UnitTests/seq.txt', str(k), str(m)]
        cmd = ["./bin/testKmerMini"] + inputs
        out = open("UnitTests/outfile.txt", "w")
        start_time = time.time()
        subprocess.run(cmd, stdout=out)
        exec_time = time.time() - start_time
        print(f"Got cpp output in {exec_time}s")
        out.close()
        parse_out("UnitTests/outfile.txt")
        os.remove("UnitTests/outfile.txt")
        os.remove("UnitTests/seq.txt")
        print("End of Minimiser")


    def test_research(self):

        def convert_to_full(superkmer:str):
            return superkmer.replace("|", "AAA")

        def parse_output(filename):
            file = open(filename, "r")
            first_line = file.readline()
            listing = first_line.split(" ")
            first_sk = convert_to_full(listing[-3])
            second_sk = convert_to_full(listing[-2])
            third_sk = convert_to_full(listing[-1]).replace('\n', '')
            all_sk = [first_sk, second_sk, third_sk]
            for line in file:
                listing = line.split(" : ")
                kmer = listing[1]
                result = int(listing[2])
                if result == -1:
                    for sk in all_sk:
                        self.assertNotIn(kmer, sk)
                else:
                    self.assertIn(kmer, all_sk[int(result)])
            file.close()

        print("======= TESTING Research =======")
        cmd = ["./bin/testResearch"]
        out = open("UnitTests/research_test.txt", "w")
        start_time = time.time()
        subprocess.run(cmd, stdout=out)
        exec_time = time.time() - start_time
        print(f"Got cpp output in {exec_time}s")
        out.close()

        parse_output("UnitTests/research_test.txt")
        os.remove("UnitTests/research_test.txt")
        print("End of Research")



if __name__ == '__main__':
    unittest.main()
