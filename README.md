(1) What do these codes do?
- Automatically processes user's Single Nucleotide Polymorphism (SNP) exome data, then compiles it with the Phase 3 of 1000 Genome data (GRCh38) data to perform 2-Dimensional (2-D) and 3-Dimensional (3-D) Principal Component Analysis (PCA).

(2) What are the required files from the user to use this code?
- VCF file containing the SNP exome data

(3) How to use these codes?
- Step 1: Upload the Python notebook (.ipynb file) to the IDE of your choice.
- Step 2: Upload your VCF file alongside all other given Python scripts (.py files) and csv files.
- Step 3: Run the Python notebook (.ipynb file).
- Step 4: Provide the necessary responses (inputs) as requested by the program.
- Step 5: Allow some time for the code to finish running (can take up to 8 hours or more depending on the hardware specifications of the user's device).

(4) What are the hardware and software requirements to run these codes?
- Python IDE (these codes were developed and tested using the standard version of Google Collaboratory with 12.67 GB RAM and 107 GB disk space)
- Python 3

(5) Why use this program?
- This program was fully developed in Python with the aim of providing a "one-stop" solution for users to extract, process and analyze their exome data alongside 1000 Genome data via PCA.
- Traditionally, bioinformaticians utilize softwares such as BCFTools and PLINK to analyze exome data from VCF files and perform PCA, respectively.
- Since this program fully runs on Python with relatively small RAM and storage requirements, users can directly execute these codes via online Python IDEs such as Google Collaboratory without having to download Jupyter, etc.

(6) What are the limitations of this program?
- This program doesn't analyze genetic mutations of the INDEL type (SNPs of INDEL type are filtered out as part of the data processing steps)
- The data compiled using the user's exome data and 1000 Genome data only undergoes one step of post-processing, namely Linkage Disequilibrium (LD) filtration.
