# ComputeCorrectKmers
#### Installation 
1. Clone the [GitHub repo](https://github.com/SaraEl-Metwally/ComputeCorrectKmers), e.g. with `git clone https://github.com/SaraEl-Metwally/ComputeCorrectKmers.git`
2. Run `make` in the repo directory for **k <= 31**  or `make k=kmersize` for **k > 31**, e.g. `make k=49`. 

#### Quick usage guide
``` ./ComputeCorrectKmers -k [kmer size] [reference genome in fasta format] [LightAssembler solid kmers file in binary format] ``` 

``` 
* [-k] kmer size                [default: 31]
* [reference kmers file in fasta format]
* [LightAssembler solid kmers file in binary format]
``` 
#### How to obtain solid kmers file from LightAssembler
1. Clone the [GitHub repo](https://github.com/SaraEl-Metwally/LightAssembler), e.g. with `git clone https://github.com/SaraEl-Metwally/LightAssembler.git`
2. Do a single-line comment `//` for the line number 945 in KmersScanning.hpp file.
3. Run `make` in the repo directory for **k <= 31**  or `make k=kmersize` for **k > 31**, e.g. `make k=49`. 
4. After running LightAssembler with your dataset, you will find the binary file of solid kmers in LightAssembler directory.
5. The name of the binary file will be [output prefix file name, i.e. default: LightAssembler].solid_kmers.

#### Outputs
The output of ComputeCorrectKmers is:
- input reference length
- number of all kmers in the reference genome 
- number of unique kmers in the reference genome
- number of all trust classified kmers by LightAssembler in Bloom filter B
- number of unique kmers in Bloom filter B
- number of incorrect kmers that are kept in Bloom filter B
- number of correct classified kmers by LightAssembler
- number of missing correct kmers by LightAssembler

#### Example 1
``` ./ComputeCorrectKmers -k 31 ```[Staphylococcus_aureus.fasta](http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/genome.fasta) [LightAssembler.solid_kmers](https://www.dropbox.com/s/jqjbcpng5u6n3mu/LightAssembler.solid_kmers?dl=0)

```
--- input reference length                          = 2903303
--- number of all kmers in the reference genome     = 2903273
--- number of unique kmers in the reference genome  = 2858856
--- number of all trust classified kmers in B       = 3185616
--- number of unique kmers in B                     = 3185616
--- number of incorrect kmers                       = 337600
--- number of correct kmers                         = 2848016
--- number of missing kmers                         = 10840
```

