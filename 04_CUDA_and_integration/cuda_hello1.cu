/* File:     cuda_hello1.cu
 * Purpose:  Implement ``hello, world'' on a gpu using CUDA.  
 *           Each thread started by the call to the kernel, 
 *           prints a message.  This version can start multiple
 *           thread blocks.
 *
 * Compile:  nvcc -o hello hello.cu 
 * Run:      ./hello <number of thread blocks> <number of threads>
 *
 * Input:    None
 * Output:   A message from each thread.
 *
 * Note:     Requires an Nvidia GPU with compute capability >= 2.0
 *
 * IPP2:     6.6 (pp. 299 and ff.)
 *
 * Sam Siewert - modified to have more friendly command line processing
 *
 */

#include <stdio.h>
#include <cuda.h>   /* Header file for CUDA */

/* Device code:  runs on GPU */
__global__ void Hello(void) {

   printf("Hello from thread %d in block %d\n", 
         threadIdx.x, blockIdx.x);
}  /* Hello */


/* Host code:  Runs on CPU */
int main(int argc, char* argv[]) {  
   int blk_ct=1;                /* Number of thread blocks */
   int th_per_blk=1;     /* Number of threads in each block */

   if(argc == 1)
   {
       printf("will use default 1 block, with 1 thread\n");
   }
   else if(argc == 2)
   {
       blk_ct = strtol(argv[1], NULL, 10);  /* Get number of blocks from command line */
   }
   else if(argc == 3)
   {
       blk_ct = strtol(argv[1], NULL, 10);  /* Get number of blocks from command line */
       th_per_blk = strtol(argv[2], NULL, 10);  /* Get number of threads per block from command line */
   }

   cudaDeviceProp prop;
   cudaGetDeviceProperties(&prop, 0); // device 0
   printf("%s, CC %d.%d\n", prop.name, prop.major, prop.minor);
   printf("maxThreadsPerBlock: %d\n", prop.maxThreadsPerBlock);
   printf("maxThreadsDim: %d x %d x %d\n",
          prop.maxThreadsDim[0],
          prop.maxThreadsDim[1],
          prop.maxThreadsDim[2]);
   printf("maxGridSize: %d x %d x %d\n",
          prop.maxGridSize[0],
          prop.maxGridSize[1],
          prop.maxGridSize[2]);
   printf("multiProcessorCount: %d\n", prop.multiProcessorCount);

   Hello <<<blk_ct, th_per_blk>>>();
            /* Start blk_ct*th_per_blk threads on GPU, */

   cudaDeviceSynchronize();       /* Wait for GPU to finish */

   return 0;
}  /* main */

