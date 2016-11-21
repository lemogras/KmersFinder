# KmersFinder

In the KmerFreq_Time project folder; There is so called lazy implementation. 
If you have enough memory for huge data files, you can use this project. 
It stores all the k-mers in RAM and makes all the operations using RAM. 
It is very fast for small files (~100-500mb) i.e. execution time is around 7-10 seconds for 100mb FASTQ files. 
However, it requires large amount of RAM for larger files. 
Even if it encodes A, T, C, G to 2-bit representation in order to store the DNA sequences,
it still requires more ram storage for files larger than 1 GB.

The program can be run via command line:
./KmersFinder_Time FileName k f

File Name : The FASTQ File
k : Lenght of k-mers
f : Most frequent f k-mers

The output of the program is list of f most frequent k-mers and their occurencies.

It requires no external libraries.
