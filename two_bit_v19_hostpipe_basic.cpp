// Owen Lucas
// Based on code from Kyle Buettner
// Two bit counting bloom filter
// Uses Algorithm from BFCounter and NEST



#include <omp.h>			// Included for timing
#include <string>			
#include <iostream>
#include <bits/stdc++.h>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include "fasta.c"

#include <sycl/ext/intel/fpga_extensions.hpp>


#include <host_pipes.hpp>
#include <CL/sycl.hpp>
using namespace cl::sycl;

#include "dpc_common.hpp"

#define MAX_SEQ_L 512

// Key parameters for Bloom filter -> These determine false positive probability (along with N)
#define M (1024*1024*512)
#define K 31
#define MAX_K_BYTES 8

#define H 8


class produce_reads1;
class parse_reads1;
class add_to_bf;

class produce_reads2;
class parse_reads2;
class query_bf;
        

struct sequence {
    std::bitset<2*MAX_SEQ_L> seq = 0;
    unsigned int             L = 0;
    bool                     valid = false;
};

struct query {
    std::bitset<MAX_SEQ_L - K + 1> found = 0;
    short seen = 0;
    unsigned int L = 0;
    unsigned int seq_num = 0;
    bool error = false;
};

struct hash_q {
    unsigned int results[H];
    unsigned int j = 0;
    unsigned int i = 0;
    unsigned int L = 0;
};

using seq_pipe = ext::intel::pipe<class some_pipe, sequence, 16>;
using seq_pipe2 = ext::intel::pipe<class other_pipe, sequence, 16>;
using res_pipe = ext::intel::pipe<class other_pipe2, unsigned int, 64>;
using res_pipe2 = ext::intel::pipe<class other_pipe3, hash_q, 16>;

//host pipe stuff

class MyPipeT;

using device2host = sycl::ext::intel::prototype::pipe<
    MyPipeT,    // An identifier for the pipe
    query,        // The type of data in the pipe
    8          // The capacity of the pipe
>;

class MyPipe2;
using host2device = sycl::ext::intel::prototype::pipe<
    MyPipe2,    // An identifier for the pipe
    sequence,        // The type of data in the pipe
    8          // The capacity of the pipe
>;

class MyPipe3;
using host2device_query = sycl::ext::intel::prototype::pipe<
    MyPipe3,    // An identifier for the pipe
    sequence,        // The type of data in the pipe
    8          // The capacity of the pipe
>;




unsigned int MurmurHash2 ( const void * key, unsigned int seed )
{
	// 'm' and 'r' are mixing constants generated offline.
	// They're not really 'magic', they just happen to work well.

	const unsigned int m = 0x5bd1e995;
	const int r = 24;
    const int length =  MAX_K_BYTES; 
	// Initialize the hash to a 'random' value

	unsigned int h = seed ^ length;

	// Mix 4 bytes at a time into the hash

	const unsigned char * data = (const unsigned char *)key;

	//while(len >= 4)
    for (int x = 4; x <= length; x = x + 4)
	{
		unsigned int k = *(unsigned int *)data;

		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h *= m; 
		h ^= k;

		data += 4;
		//len -= 4;
	}
	
	// Handle the last few bytes of the input array
    
    int rem = length & 3;
    
    if (rem != 0){

        switch(rem)
        {
        case 3: h ^= data[2] << 16;
        case 2: h ^= data[1] << 8;
        case 1: h ^= data[0];
                h *= m;
        };
    }

	// Do a few final mixes of the hash to ensure the last few
	// bytes are well-incorporated.

	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;

	return h;
}


 

int main (int argc, char *argv[]) {
    
    
    #if defined(FPGA_EMULATOR)
        ext::intel::fpga_emulator_selector device_selector;
    #else
        ext::intel::fpga_selector device_selector;
    #endif
    
    try{
        FILE *database_file;
	
        FASTAFILE *ffp;
        
        char *seq;
        char *name;
        int L;

        char* file_name = "./data/sra_data.fasta";

        ffp = OpenFASTA(file_name);
        long long unsigned int num_database_char = 0;

        char c;

        char kmer[K+1];
        kmer[K] = '\0';
        long long unsigned int reverse_mask = pow(2, (2 * K));
        reverse_mask = reverse_mask - 1;
        long long unsigned int one = 1;
        long long unsigned int check_mask = pow(2, (2 * K -1)); //check if the top bit is 1
        
    
        int total_valid_reads = 0;
        unsigned long long int total_nucleotides = 0;
        unsigned long long int total_kmers = 0;

        const int max_reads = 100;

        
        sequence* all_reads = (sequence*)malloc(sizeof(sequence) * max_reads);
        
        query* all_qs = (query*)malloc(sizeof(query) * max_reads);
        
        query* all_qs_HP = (query*)malloc(sizeof(query) * max_reads);
        
        short* kmer_len = (short*)malloc(sizeof(short) * max_reads);
        
        
        //read the sequecing reads from the file
         while (ReadFASTA(ffp, &seq, &name, &L) && total_valid_reads < max_reads) {

            if (L >= K){
                
                if(L > MAX_SEQ_L){
                    L = MAX_SEQ_L;
                }
                
                
                std::bitset<2*MAX_SEQ_L> seq_bitset = 0;
                std::bitset<2*MAX_SEQ_L> char_as_bits = 0;

                for (int i = 0; i < L; i++){
                    char_as_bits = 0; //any unknown char gets input as an A, should figure out what to actually do
                    if (seq[i] == 'a' || seq[i] == 'A'){
                        char_as_bits = 0;
                    }
                    else if(seq[i] == 'c' || seq[i] == 'C'){
                        char_as_bits = 1;
                    }
                    else if(seq[i] == 'g' || seq[i] == 'G'){
                        char_as_bits = 2;
                    }
                    else if(seq[i] == 't' || seq[i] == 'T'){
                        char_as_bits = 3;
                    }
                    seq_bitset <<= 2;
                    seq_bitset |= char_as_bits;
                }
                
                all_reads[total_valid_reads].seq = seq_bitset;
                all_reads[total_valid_reads].L = L;
                all_reads[total_valid_reads].valid = true;
                
                kmer_len[total_valid_reads] = (L- K + 1);
                
                total_valid_reads++;
                total_nucleotides += L;
                total_kmers += (L- K + 1);
            }
            free(seq);
            free(name);
        }
        
        CloseFASTA(ffp);

        
        // Create a queue bound to the chosen device.
        // If the device is unavailable, a SYCL runtime exception is thrown.
        printf("Going to attempt to create the queue\n");
        sycl::queue q(device_selector, dpc_common::exception_handler);
        printf("Queue has been made\n");

        // Print out the device information.
        std::cout << "Running on device: "
                  << q.get_device().get_info<info::device::name>() << "\n";

        bool compose_finished = false;


        // Phase 1.a - Read sequences from host pipe, write to phase1.b
        auto e = q.submit([&](handler &h) {
            h.single_task<produce_reads1>([=]() {
                sequence seq;
                for (int j = 0; j < total_valid_reads; j++){
                    seq = host2device::read();
                    bool success = false;
                    do{
                        seq_pipe::write(seq, success);
                    }while(!success);
                }
            });
        });

                
        // Phase 1.b - read from phase 1.a, split up data an write smaller chunks to phase 1.c
         auto e2 = q.submit([&](handler& h) {
            h.single_task<parse_reads1>([=]()[[intel::kernel_args_restrict]]{
                for (int j = 0; j < total_valid_reads; j++){
                    sequence dev_seq = seq_pipe::read();
                    for(int i = 0; i < dev_seq.L - K + 1; i++){
                        std::bitset<2*K> kmer_bits = 0;
                        std::bitset<MAX_K_BYTES * 8> kmer_bytes = 0;
                        for(int n = 0; n < H; n++){
                            unsigned int result = MurmurHash2(&kmer_bytes, n) & (unsigned long long int)(M-1);
                            //write the result to the pipe
                            res_pipe::write(result);
                        }
                    }
                }
            });
        });
        
        //Phase 1.c - read from phase 1.b
        auto acc_bf = q.submit([&](handler& h) {
            h.single_task<add_to_bf>([=]()[[intel::kernel_args_restrict]]{
                for (int j = 0; j < total_kmers*H; j++){
                    unsigned int result = res_pipe::read();
                }
            });
        });
        
        printf("\nGoing to start writing sequences to the device\n");
        
        for (int j = 0; j < total_valid_reads; j++){
            host2device::write(q, all_reads[j]);
        }
        printf("\nDone writing sequences to the device\n");
        

        //Begin phase 2

        // Phase 2.1 - Read from second host pipe, write to 2.b
        auto e3 = q.submit([&](handler &h) {
            h.depends_on(acc_bf);
            h.single_task<produce_reads2>([=]() {
                sequence seq;
                for (int j = 0; j < total_valid_reads; j++){
                    seq = host2device_query::read();
                    bool success = false;
                    do{
                        seq_pipe2::write(seq, success);
                    }while(!success);
                }

            });
          });

        // Phase 2.b - read from 2.a, split up data into smaller chunks and write to 2.c
        q.submit([&](handler& h) {
                h.single_task<parse_reads2>([=]()[[intel::kernel_args_restrict]]{
                    for (int j = 0; j < total_valid_reads; j++){
                        sequence dev_seq = seq_pipe2::read();
                        for(int i = 0; i < dev_seq.L - K + 1; i++){
                            hash_q curr_res;
                            curr_res.i = i;
                            curr_res.j = j;
                            curr_res.L = dev_seq.L;
                            res_pipe2::write(curr_res);
                        }
                    }
                });
            });
        
        // Phase 2.c - read small chunks from 2.b, write results back to thrird hostpipe
        auto acc_bf2 = q.submit([&](handler& h) {
            h.single_task<query_bf>([=]()[[intel::kernel_args_restrict]]{
                query query_buf[16]; 
                for (int j = 0; j < total_kmers; j++){
                    hash_q read_result = res_pipe2::read();
                    std::bitset<H> found;
                    found.reset();
                    
                    unsigned int idx = read_result.j & 15;
                    if (query_buf[idx].seen == 0){//first time writing to the buffer after clear
                        query_buf[idx].seq_num = read_result.j;
                    }
                    else if(query_buf[idx].seq_num != read_result.j){ //circle buffer has been overflowed
                        query_buf[idx].error = true;
                    }
                    
                    query_buf[idx].seen++;
                    
                    if (query_buf[idx].seen == (read_result.L - K + 1)){
                        device2host::write(query_buf[idx]);
                        query_buf[idx].seen = 0;
                        query_buf[idx].error = false;
                        query_buf[idx].found = 0;
                    }
                }
            });
        });
        
        
        printf("\nStarting to read from the host pipe\n");
        
        //need to handle both reading and writing from the board at the same time
        bool write_success = false;
        bool read_success = false;
        int write_count = 0;
        int read_count = 0;
        
        
        while (write_count < total_valid_reads || read_count < total_valid_reads){
            if(write_count < total_valid_reads){//try to write if there are still more sequences to send
                host2device_query::write(q, all_reads[write_count], write_success);
                if(write_success){
                    write_count++;
                    write_success = false;
                }
            }
            if(read_count < total_valid_reads){ //try to read if there are still more results to get back
                query dev_q = device2host::read(q, read_success);
                if(read_success){
                    all_qs_HP[dev_q.seq_num] = dev_q;
                    
                    read_count++;
                    read_success = false;
                }
            }
        }
        
        
        printf("\nDone reading from the host pipe\n");

        free(all_reads);
        free(all_qs);
        free(kmer_len);
       
        printf("Program complete\n");
        return 0;
    }
    catch (sycl::exception const& e) {
    // Catches exceptions in the host code
    std::cerr << "Caught a SYCL host exception:\n" << e.what() << "\n";

    // Most likely the runtime couldn't find FPGA hardware!
    if (e.code().value() == CL_DEVICE_NOT_FOUND) {
      std::cerr << "If you are targeting an FPGA, please ensure that your "
                   "system has a correctly configured FPGA board.\n";
      std::cerr << "Run sys_check in the oneAPI root directory to verify.\n";
      std::cerr << "If you are targeting the FPGA emulator, compile with "
                   "-DFPGA_EMULATOR.\n";
    }
    std::terminate();
  }

}


