# sk7
This repository is a library made to manipulate sorted sets of super-k-mers, and is made to be compatible with Kff files.
The code is not guaranty to be bug-free.

## Install
```shell
git clone https://github.com/yoann-dufresne/sk7.git
git submodule init
git submodule update
```

Be sure to use the branch debug_compact / compact of the kff-tools submodule,
to ensure minimizer compatibility.

## How to use the library

All the library content is in a namespace named sk7 (except Kff_scanner).

### Input

Here is an example (from apps/main.cpp) to read the kmers from a kff file:

```C++
Kff_scanner* scanner = new Kff_scanner(input_name, m, bucketed, sorted);
BucketMap* read = scanner->readAll();
delete scanner;
```

where:
- input_name is the name of the kff file.
- m is the wanted size of the minimizers.
- bucketed is to be set to true if the kff file went through `kff-tools bucket`.
- sorted is to be set to true if the kff file went through `kff-tools compact -s`.

This will create a BucketMap object that can search Kmer objects with the find method.

The constructor of Kff_scanner will initiate the library with the values of k and m found and/or provided.

### Manual initialisation

The following function can be used to set the values of k and m in the classes:

```C++
    void initLib(int _k, int _m, bool (*_infKmer) (const Kmer &, const Kmer &) = &infId, bool (*_equalKmer) (const Kmer &, const Kmer &) = &equalId);
```

where:
- _k is the length of the k-mers.
- _m is the wanted size of the minimizers.
- _infKmer (optional) add the possibility to define a custom comparison function for Kmers.
- _equalKmer (optional) add the possibility to define a custom equality function for Kmers.

### Classes description

This section will consist of a short presentation of the classes used in this library, the code is commented for more
details.

#### Kmer

The content of a Kmer is represented by an uint128_t so Kmers cannot represent every size possible.
A Kmer of length 0 is meant to represent a non-existent Kmer (necessary for an aspect of the binary search in a Bucket_).

#### SuperKmer

The content of a SuperKmer (see [here](https://hal.archives-ouvertes.fr/hal-02435086/document) for definition), is 
represented as a vector of uint8_t. The nucleotide are sorted in interleaved order, here is an example:

Let's take k = 5 and m = 3:
GCAAAG -> minimizer = AAA so we store <span style="color:green">G</span><span style="color:red">C</span>
|<span style="color:blue">G</span> in this order <span style="color:blue">G</span><span style="color:red">C</span>_
<span style="color:green">G</span>

The first `ceil(log2(k - m + 1))` bits of the vector are used to store the prefix and suffix length (in nucleotides), 
the rest stores the data. You can use the getter functions to read directly those values (the getValue return an uint64_t,
be careful to the size of the data).

The print function prints the nucleotides in interleave order too.

#### Minimizer

A class that contains a m-mer stored as an uint64_t. The constructor takes a function as a parameter that must return a 
struct HashPos with the value of the found minimizer and its position (the default function alpha is available in 
exampleHash.cpp), allowing to compare with the reverse complement.

#### Bucket_

A Bucket_ is a set (represented as vector of SuperKmer). It's support binary search of a Kmer.
The addToList method doesn't guaranty the order, so to add a Kmer use AddKmer or AddSuperKmer.
It supports intersection, union and symmetrical difference that break the compaction in SuperKmer.
The chainedUnion method is a union that keep compaction but is still experimental.
The print method allows to print every SuperKmer in the bucket in the order they're stored.

#### BucketMap

A BucketMap index every Bucket_ contained in its map attribute. It can search for a Kmer or add a Kmer in the right
Bucket_.