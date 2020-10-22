# The Simple User Manual for TBpipeline

1. Get the repository from the GitHub and change your directory to the local repository.

   ```bash
   # Clone the Assignment1 Repository
   $ git clone https://github.com/B176888-2020/Assignment1.git
   
   # Change the directory to your local Assignment1 directory
   $ cd ./Assignment1
   ```

2. Unzip the B176888-2020.Assignment1.zip with  `[password]` in the B176888-2020.Assignment1.pdf. The bash scirpt, `TBpipeline.sh` and git logfile, `all_the_things_I_did` will be revealed.

   ```bash
   $ unzip -P [password] B176888-2020.Assignment1.zip
   ```

3. (Optional)By default, the `TBpipeline.sh` has already be set as an executable program. If TBpipeline.sh is not executable, then try the following command to change its mode. 

   ```bash
   $ chmod +x ./TBpipeline.sh
   ```

4. For this assignment, you can run these commands:

   **Attention**: The directories arguments should include the final "/" as it is convenient for auto-completion and they are directories containing targeted files rather than the path to the files themselves.

   1. Run the whole pipeline without a pause.

      ```bash
      $ ./TBpipeline.sh -f /localdisk/data/BPSM/Assignment1/fastq/ -r /localdisk/data/BPSM/Assignment1/Tbb_genome/ -b /localdisk/data/BPSM/Assignment1/ -o ./output/ -y
      ```

   2. Run the pipeline with a a breakpoint to decide whether to continue downstream processes based on the FastQC results' summary

      ```bash
      $ ./TBpipeline.sh -f /localdisk/data/BPSM/Assignment1/fastq/ -r /localdisk/data/BPSM/Assignment1/Tbb_genome/ -b /localdisk/data/BPSM/Assignment1/ -o ./output/
      ```

5. More general usage and parameters chart will be listed below

   ```bash
   $ ./TBpipeline.sh -f [fqfile_and_fq.gz_dir] -r [ref_genome_dir] -b [bedfile_dir] -o [output_dir] -y
   ```

   | Parameters | Reasons/Usage                                                |
   | ---------- | ------------------------------------------------------------ |
   | -f         | The directory of fqfile and  fq.gz files                     |
   | -r         | The directory of the reference  genome data                  |
   | -o         | The directory to stored your  outputs                        |
   | -y         | The confirmation to perform  downstream analysis after FastQC |

   

