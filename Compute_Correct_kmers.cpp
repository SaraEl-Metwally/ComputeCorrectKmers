//Usage: ./ComputeCorrectKmers -k <kmer_size> <genome file in fasta > <solid kmers file resulted from LightAssembler>
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include <cstdlib>
#include <map>
#include <vector>
#include "BinaryStore.hpp"
#include <getopt.h>//getopt_long
#include <sys/stat.h>//stat
#include <sstream>
#ifdef largeintlib

#include "LargeInt.hpp"

typedef LargeInt<kmer_precision> kmercode_length;

#else

#if (! defined kmercode_length) || (! defined _LP64)

typedef uint64_t kmercode_length;

#endif
#endif
using namespace std;

struct program_options

{

    int k;

    std::vector<std::string> read_files;

    program_options():k(31){}

};

bool check_options(program_options &opt)

{

    bool success=true;

    std::cerr<<std::endl;

    //************************** kmer size**************************************************
   
    if(opt.k>0)

    {

       if(opt.k>((int)sizeof(kmercode_length)*4))

          {

          

           std::cerr<<"--- maximum support kmer size for this compiled version : "<<(sizeof(kmercode_length)*4)<<std::endl;

           std::cerr<<"--- use make k="<<opt.k<<std::endl;

           success=false;

         

          }

       else if(opt.k%2==0)

          {opt.k--;std::cout<<"--- to avoid palindromes, kmer size must be odd, suggested kmer size "<<opt.k<<std::endl;}

    }

    else

    {std::cerr<<"--- invalid value for kmer size: "<<opt.k<<std::endl;success=false;}


    if(opt.read_files.size()==0)

    {std::cerr<<"--- no read files specified as inputs"<<std::endl;success=false;}

    else

    {

        struct stat stat_file_info;

        int int_stat;

        std::vector<std::string>::const_iterator it;

        for(it=opt.read_files.begin();it != opt.read_files.end();++it)

        {

            int_stat=stat(it->c_str(),&stat_file_info);

            if(int_stat != 0)

            {

                std::cerr<<"--- error: file not found "<<*it<<std::endl;
    
                success=false; 
            }
        }


      }

   return success;



}




int nt2int(char nt)

{

    

    if(nt=='A'||nt=='a')

        return 0;

    if(nt=='C'||nt=='c')

        return 1;

    if(nt=='G'||nt=='g')

        return 2;

    if(nt=='T'||nt=='t')

        return 3;

}





kmercode_length get_reverse_complement(kmercode_length code, int kmer_size)
{
    int i ;
    kmercode_length rev_comp =  (static_cast<kmercode_length>(0)) ;
    for ( i = 0 ; i < kmer_size ; ++i )
    {
        kmercode_length tmp = ( code >> ( 2 * i ) ) & (static_cast<kmercode_length>(3)) ;
        rev_comp = ( rev_comp << 2 ) | ( (static_cast<kmercode_length>(3)) - tmp ) ;
    }


    return rev_comp;
}

kmercode_length get_canonical_kmer_code(kmercode_length code, int kmer_size)
{
    int i ;
    kmercode_length rev_com = get_reverse_complement(code, kmer_size);
    return rev_com < code ? rev_com : code ;
}




void print_usage()

{

    std::cerr <<std::endl;

    std::cerr <<" ******************** <<<  Accuracy of LightAssembler  >>> *********************** " <<std::endl<<std::endl ;

    std::cerr <<" Compute the number of correct classified kmers compared to the reference genome."<<std::endl<<std::endl ;

    std::cerr <<" Usage: ./ComputeCorrectKmers -k < kmer size > < reference genome > < LightAssembler solid kmers file > "<<std::endl;

    std::cerr <<std::endl<<

    "--- [-k] kmer size [default: 31] "<<std::endl<<

    "--- < reference genome in fasta format >  "<<std::endl<<

    "--- < LightAssembler solid kmers file in binary format > "<<std::endl<<std::endl;
    

    std::cerr <<" Typical ComputeCorrectKmers Command Line :"<<std::endl<<std::endl;

    std::cerr <<" ./ComputeCorrectKmers -k 31 reference_genome.fasta LightAssembler.solid_kmers "<<std::endl;

    std::cerr <<std::endl;



}

void parse_options(int argc,char **argv,program_options &opt)

{

   const char* opt_string =":k:";

   static struct option long_options[]=
   {
      {"kmer size",required_argument,0,'k'},
      {NULL,no_argument,NULL,0}

   };
   int option_index=0;

   int ch;

   while((ch=getopt_long(argc,argv,opt_string,long_options,&option_index))!=-1)
   {

        if(ch=='k')
         {     
             std::istringstream argument(optarg);
             if(!(argument>>opt.k)||!(argument.eof()))
             {
                std::cerr<<std::endl;
                std::cerr<<"--- invalid argument of -k "<<optarg<<std::endl;
                print_usage();
                exit(1);
             }

         }
        else if(ch=='?')
        {

             std::cerr<<std::endl;
             if (isprint (optopt))
             std::cerr<<"--- invalid option -"<< static_cast<char>(optopt) <<std::endl;
             else
             std::cerr<<"--- invalid option "<< argv[optind-1]<<std::endl;
             print_usage();
             exit(1);
        }
        else if(ch==':')
        {

             std::cerr<<std::endl;
             std::cerr<<"--- missing argument of -"<< static_cast<char>(optopt)<<std::endl;
             print_usage();
             exit(1);

        }

 }

   for(int i=optind;i<argc;++i)

   {opt.read_files.push_back(argv[i]);}

}



int main(int argc,char** argv)
{
 
   program_options opt;

   parse_options(argc,argv,opt);

   if(!check_options(opt))

    {

        print_usage();

        exit(1);

    }

 int kmer_size=opt.k;
 kmercode_length kmercode=0,kmercode_can=0;
 string file_ref=opt.read_files[0];
 string file_kmers=opt.read_files[1];
 ifstream file_kmers_Ref(file_ref.c_str(),std::ios::in);
 string solid_kmers_file(file_kmers); 
 binary_store solid_kmers(solid_kmers_file,sizeof(kmercode),false);
 std::map<kmercode_length, unsigned int> hash_table_ref;
 std::map<kmercode_length, unsigned int> hash_table_B;
 uint64_t nb_ref_kmers=0;
 uint64_t nb_ref_kmers_all=0;
 uint64_t nb_kmers_B=0;
 uint64_t nb_kmers_B_all=0;
 uint64_t nb_correct_kmers=0;
 uint64_t nb_missing_kmers=0;
 uint64_t nb_incorrect_kmers=0;
 kmercode_length kmer_mask=((static_cast<kmercode_length>(1))<<(kmer_size*2))-1;
 if (file_kmers_Ref.is_open())
    {
       string read = "";
       getline(file_kmers_Ref, read);
       read ="";
       string ref="";
       while (file_kmers_Ref >> read)
       {
            ref= ref+read;

       }

       file_kmers_Ref.close();
       int i;
       std::cout<<"--- input reference length     = "<<ref.length()<<std::endl;
       kmercode=0,kmercode_can=0;
  
       for(i=0; i<kmer_size; ++i)
       {
          kmercode=kmercode*4+nt2int(read[i]);
       }
          kmercode_can=get_canonical_kmer_code(kmercode, kmer_size);
          hash_table_ref[kmercode_can]=1;
          nb_ref_kmers_all++;
       for(i=1; i<ref.length()-kmer_size+1; ++i)
       {
          kmercode=(kmercode*4+nt2int(ref[i+(kmer_size-1)])& kmer_mask);
          kmercode_can=get_canonical_kmer_code(kmercode, kmer_size);
          nb_ref_kmers_all++;
          if (hash_table_ref.find(kmercode_can) == hash_table_ref.end() )
            {

                hash_table_ref[kmercode_can]=1;

                nb_ref_kmers++; //number of unique reference kmers without repetition.

            }

          else

                hash_table_ref[kmercode_can]=hash_table_ref[kmercode_can]+1;         
        }


    }
 else
    {
        std::cerr << "--- can't open the input Genome reference file " << std::endl;
        exit(1);

    }


    std::cout<<"--- number of all kmers in the reference genome     = "<<nb_ref_kmers_all<<std::endl;
    std::cout<<"--- number of unique kmers in the reference genome     = "<<nb_ref_kmers<<std::endl;
    solid_kmers.rewind_all();
    kmercode=0;
    while(solid_kmers.read_element(&kmercode))
    {

            nb_kmers_B_all++;
            if (hash_table_B.find(kmercode) == hash_table_B.end() )

            {

                hash_table_B[kmercode]=1;

                nb_kmers_B++; //number of unique kmers in B without repetition.

            }



            else

            hash_table_B[kmercode]=hash_table_B[kmercode]+1;
        

    }


    std::cout<<"--- number of all trust classified kmers in B = "<<nb_kmers_B_all<<std::endl;
    std::cout<<"--- number of unique kmers in B     = "<<nb_kmers_B<<std::endl;
    std::map<kmercode_length, unsigned int>::iterator itr;
    for(itr = hash_table_B.begin(); itr != hash_table_B.end(); itr++) 
    {
          if (hash_table_ref.find(itr->first) == hash_table_ref.end() )
            {

               

                nb_incorrect_kmers++; //number of incorrect kmers without repetition.

            }

          else

                nb_correct_kmers++;




    }


   std::cout<<"--- number of incorrect kmers     = "<<nb_incorrect_kmers<<std::endl;
   std::cout<<"--- number of correct kmers       = "<<nb_correct_kmers<<std::endl;
   nb_missing_kmers= nb_ref_kmers - nb_correct_kmers;
   std::cout<<"--- number of missing kmers       = "<<nb_missing_kmers<<std::endl;


 return 0;
}
 
















