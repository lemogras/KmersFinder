//
//  main.cpp
//  KmersFreq_Time
//
//  Created by Baki Er on 20.11.2016.
//  Copyright © 2016 Baki Er. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;


int main(int argc, const char * argv[]) {
    
    string line;                        // For reading each line of FastQFile
    string path_file;                   //Path of FastQFile
    stringstream kValue,fValue,path;    // For command line arguments
    int k,f;
    
    std::ifstream FastQFile;            // Input FastQ File
    ifstream SeqBinTxt;                 // These three are for reading and writing batch of k-mers from/to hdd.
    ofstream SeqBinTxtSorted;
    ofstream SeqBinTxtSortedFinal;
    
    
    vector<unsigned long long int> kMerSeqBin;                          // A vector containing k-mers' binary encoding
    vector<std::pair<long long unsigned int,int>> KmersFreq;            // For counting kmers and corresponding counts
    vector<std::pair<int,long long unsigned int>> FrequentKmersResult; // Desired "f" most frequent kmers and counts
    
    int N;                      // Length of each DNA Sequence
    
    unsigned long long binrep; // These are for encoding A,T,C,G to binary
    unsigned long long tempbin;
    char tempNuc;
    string genomeFragment;
    
    
    //Some dummy variables used through main function
    int count = 1;
    unsigned long long int temp=0;
    bool initialFlag;
    
    // Command Line argument assignments
    kValue << argv[2];
    fValue << argv[3];
    path << argv[1];
    kValue >> k;
    fValue >> f;
    path >> path_file;
    
    FastQFile.open(path_file);              //Opening Fast Q file

    //READING DNA SEQUENCES BATCH BY BATCH FROM THE FASTQFILE
    
    while(!FastQFile.eof())
    {
        getline(FastQFile,line);
        getline(FastQFile,line);
        N = int(line.length());
        for (int m=0; m<N-k+1 ;m++)
        {
            genomeFragment =line.substr(m,k);
            binrep = 0b00;
            for (int t=0; t<k; t++)
            {
                tempNuc = genomeFragment[t];
                tempbin =  tempNuc & 0b00000110;
                tempbin = (tempbin >> 1);
                binrep = (binrep<<2)|tempbin;
            }
            kMerSeqBin.push_back(binrep);
        }
        getline(FastQFile,line);
        getline(FastQFile,line);
    }
    
    FastQFile.close();
    
    //SORTING ALL KMERS
    
    std::sort(kMerSeqBin.begin(),kMerSeqBin.end());
    
    
    //COUNTING OCCURENCIES OF A KMERS AND STORING WITH FREQUENCIES

    
    count = 0;
    initialFlag = true;
    for (const auto &y: kMerSeqBin)
    {
        if(initialFlag) temp = y; else;
        if(temp == y)
        {
            count++;
            initialFlag = false;
        }
        else
        {
            FrequentKmersResult.push_back(std::pair<int,long long unsigned int>(count,temp));
            count = 1;
            temp = y;
        }
    }
    
    kMerSeqBin.clear();
    
    //SORTING ACCORDING TO THEIR FREQ. AND PRINTING MOST FREQUENT ONES
    
    std::sort(FrequentKmersResult.begin(), FrequentKmersResult.end());
    
    //DECODING BEFORE PRINTING
    for(int i = int(FrequentKmersResult.size())-1; i > int(FrequentKmersResult.size())-1-f; i--)
    {
        
        std::bitset<128> seqbin(FrequentKmersResult[i].second);
        std::string FreqKmer;
        
        for (int b=k*2-1; b>-1; b=b-2)
        {
            if(seqbin[b]==0)
            {
                if (seqbin[b-1]==0)
                {
                    FreqKmer +='A';
                }
                else FreqKmer +='C';
            }
            else
            {
                if (seqbin[b-1]==0)
                {
                    FreqKmer +='T';
                }
                else FreqKmer +='G';
                
            }
        }
        std::cout << FreqKmer << " frequency =" << FrequentKmersResult[i].first << std::endl;
    }
    return 0;
}

