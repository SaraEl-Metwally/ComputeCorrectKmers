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
