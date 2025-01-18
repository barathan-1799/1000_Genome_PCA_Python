(1) What's the use of this program?
- Automatically extracts and processes user's Single Nucleotide Polymorphism (SNP) exome data, then compiles it with the Phase 3 of 1000 Genome data (GRCh38) data based on matching chromosomal loci (POS) of the SNPs to generate 2-Dimensional (2-D) and 3-Dimensional (3-D) Principal Component Analysis (PCA) plots.

(2) What are the required files from the user to use this program?
- VCF file containing the SNP exome data

(3) How to use this program?
- Step 1: Upload the Python notebook (.ipynb file) to the IDE of your choice.
- Step 2: Upload the VCF file containing your exome data alongside all other given Python scripts (.py files) and csv files.
- Step 3: Run the Python notebook (.ipynb file).
- Step 4: Provide the necessary responses (inputs) as requested by the program.
- Step 5: Allow some time for the code execution (can take up to 8 hours or more depending on the hardware specifications of the user's device).

(4) What are the hardware and software requirements to run this program?
- Python IDE (these codes were developed, tested and executed using the standard version of Google Collaboratory with 12.7 GB RAM and 107.7 GB disk space)
- Python 3

(5) Why use this program?
- This program was fully developed in Python with the aim of providing a "one-stop" solution for users to extract, process, compile, and perform PCA on their own exome data alongside the Phase 3 of 1000 Genome data.
- Over 40,000 SNPs can be analyzed by this program using hardware requirements that do not exceed 12.7 GB RAM and 107.7 GB disk space.
- Traditionally, bioinformaticians utilize softwares such as BCFTools and PLINK to analyze exome data from VCF files and perform PCA, respectively.
- Since this program fully runs on Python with relatively small RAM and storage requirements, users can directly execute these codes via online Python IDEs such as Google Collaboratory without having to download Jupyter, etc.
- This program directly extracts the Phase 3 of 1000 Genome data VCF files available online (https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes) without requiring the user to individually download and unzip these large-sized files on their respective devices (the zipped VCF file for each chromosome can be as large as 1.1 GB).

(6) What are the limitations of this program?
- This program doesn't analyze genetic mutations of the INDEL type (SNPs of INDEL type are filtered out as part of the data processing step).
- The data compiled using the user's exome data and 1000 Genome data only undergoes one step of post-processing, namely filtration based on a user-specified Linkage Disequilibrium (LD) threshold value.
