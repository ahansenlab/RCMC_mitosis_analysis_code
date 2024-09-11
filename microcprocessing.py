#Aligning and processing a single fastq file through a pipeline similar to distiller but instead using bowtie2

from sys import exit
import subprocess as sp
import argparse
import multiprocessing
import uuid

parser = argparse.ArgumentParser(description = "run bwa/bowtie2 and pairtools on fastq files to produce pairsam files")
bamopts = parser.add_mutually_exclusive_group()
parser.add_argument("--file_1", "-1", help = "first demuxed fastq of paired end reads - required", nargs = "*")
parser.add_argument("--file_2", "-2", help = "second demuxed fastq of paired end reads - required", nargs = "*")
parser.add_argument("--genome", "-g", help = "genome to align to - mouse or human - required")
parser.add_argument("--genometype", "-y", help = "genome type - use if your genome is a modified version of a standard genome - should be one of hg19, hg38, mm9, mm10, or mm39")
parser.add_argument("--threads", "-t", help = "number of threads to use for alignment - default is 1", default = "1")
parser.add_argument("--resolutions", "-r", help = "list of resolutions to output in decreasing order - all resolutions must be a multiple of the smallest resolution - default: 10000000 5000000 2500000 1000000 500000 250000 100000 50000 25000 10000 5000 2000 1000", nargs = "*")
parser.add_argument("--out", "-o", help = "name for output files - defaults to name of first file")
parser.add_argument("--outdir", help = "a directory to store output files - default is current directory", default = "./")
parser.add_argument("--aligner", "-a", help = "choose the aligner to use - options are bwa or bowtie2 - bowtie2 works for short (~50 bp) reads, but bwa is needed for long reads - required")
parser.add_argument("--walks", "-w", help = "choose walks setting - options are [mask] (discard these reads), [***5any (report 5' most alignment on each end),***NOT CURRENTLY AVAILABLE] [all]/[parse2] (in fragment A-B-C report A-B, B-C (not A-C)), [expand] (in fragment A-B-C report A-B, B-C, A-C) - required")
parser.add_argument("--filter", "-f", help = "Set the filter option for reads - the best choice depends on which aligner you are using: bowtie2 default is reads with MAPQ >= 2, bwa default is >= 30. The specified value can be a number for a minimum acceptable MAPQ score (eg 2, 30, 40), or it can be XS for bowtie2 (removes any read with a secondary mapping site, sets MAPQ minimum to 1)")
parser.add_argument("--captureregion", "-c", help = "[THIS DOESN'T DO ANYTHING YET :)] if analysing data from RCMC supply regions of interest here - this will also turn on RCMC mode and produce pairs, cools, and mcools which only contain the reads of interest in addition to the full files - space separate regions written like this: chrA:XXXXXXXX-XXXXXXXX chrB:XXXXXXXX-XXXXXXXX", nargs = "*")
bamopts.add_argument("--alignonly", "-b", help = "[ONLY BOWTIE2 FOR NOW] only run bowtie2/bwa and make bams - can be useful for post initial analysis QC steps", action = "store_true")
bamopts.add_argument("--keepbams", "-k", help = "keep bam files while doing a normal full analysis", action = "store_true")
args = parser.parse_args()

file1 = args.file_1
file2 = args.file_2
genome = args.genome
gentype = args.genometype
threads = args.threads
outname = args.out
outdir = args.outdir
reslist = args.resolutions
alignonly = args.alignonly
keepbams = args.keepbams
captureregion = args.captureregion
aligner = args.aligner
walks = args.walks
filter = args.filter




#Check requirements are fulfilled:
## TODO: update this to account for the new aligner
condapacks = sp.run("conda list".split(), capture_output=True)
condapacksstr = str(condapacks.stdout)
if alignonly:
    if "bowtie2" not in condapacksstr or "samtools" not in condapacksstr or "sambamba" not in condapacksstr:
        print("Please make sure bowtie2, samtools and sambamba are installed in your current conda environment (check conda list)")
        exit()
elif not alignonly:
    if "bowtie2" not in condapacksstr or "pairtools" not in condapacksstr or "cooler" not in condapacksstr or "pairix" not in condapacksstr:
        print("Please make sure bowtie2, pairtools, pairix and cooler are installed in your current conda environment (check with 'conda list')")
        exit()

#Check that outdir ends with a /, add one if it doesn't
if args.outdir is not None and not outdir.endswith("/"):
    outdir = outdir + "/"

#Check a genome was specified
if args.genome is None:
    print("Genome not specified - check help for formatting")
    parser.print_usage()
    exit()

#Check an aligner was specified
if args.aligner != "bowtie2" and args.aligner != "bwa":
    print("Aligner not specified - check help for formatting")
    parser.print_usage()
    exit()
elif args.aligner == "bowtie2":
    if filter is None:
        filterexp = "2"
        optbtfilter = ""
    elif filter == "XS":
        optbtfilter = "| grep -v XS: - "
        filterexp = "1"
    else:
        filterexp = filter
        optbtfilter = ""
elif args.aligner == "bwa":
    if filter is None:
        filterexp = "30"
        optbtfilter = ""
    elif filter == "XS":
        print("Filter option XS only works with bowtie2!")
        exit()
    else:
        filterexp = filter
        optbtfilter = ""

#Check if file1 and file2 are single files or multiple
if file1 is None or file2 is None:
    print("Input files not specified - check help for formatting")
    parser.print_usage()
    exit()
elif len(file1) > 1 and type(file1) == list and type(file2) == list and len(file1) == len(file2):
    multifile = 1
    if args.out is None:
        outlist = [fname + "_" + genome for fname in file1]
        outname = outlist[0]
    else:
        outlist = list()
        for i in range(len(file1)):
            outlist.append(outname + "_" + str(i + 1))
    #Input to pair merging step needs all of the outputs together
    pairnamelist = [outdir + oname + ".pairs" for oname in outlist]
    pairnamest = " ".join(pairnamelist)
    bamlist = [outdir + oname + ".sorted.bam" for oname in outlist]
    bamst = " ".join(bamlist)
    bailist = [outdir + oname + ".sorted.bam.bai" for oname in outlist]
    baist = " ".join(bailist)
elif len(file1) == 1 and len(file2) == 1:
    multifile = 0
    #If nargs = *, always makes a list, even if only one element
    file1 = "".join(file1)
    file2 = "".join(file2)
    if args.out is None:
        outname = file1 + "_" + genome
    pairnamest = outdir + outname + ".pairs"
else:
    print("Mismatch in number of input files, check arguments")
    exit()

#Check that a sensible number of threads has been requested - more protections here are possible - at the moment users are trusted to be sensible
cpucount = multiprocessing.cpu_count()
if args.threads is None:
    print("Defaulting to one thread")
    threads = 1
elif int(args.threads) >= cpucount:
    print("Too many threads requested, resetting to default")
    threads = 1

#Check that the user has entered a valid genome to align to
if args.genometype is None:
    gentype = args.genome

if gentype == "mm10" or gentype == "mm39" or gentype == "mm9":
    toprint = "Aligning to mouse genome {}".format(gentype)
    print(toprint)
elif gentype == "hg19" or gentype == "hg38":
    toprint = "Aligning to human genome {}".format(gentype)
    print(toprint)
else:
    if gentype == genome:
        print("Genome option not recognised or not entered. Please use mm9/10/39 or hg19/38 or ask Miles to change the script to accommodate your new organism/genome. If you are using a modified version of base genome, use the -g option to indicate the base genome name.")
        exit()
    else: #If they're using a modified genome, make sure the base genome exists so that the files are redirected properly
        print("Genome/base genome option not recognised. Please use mm9/10/39 or hg19/38 or ask Miles to change the script to accommodate your new organism/genome.")
        exit()

#Set up resolutions as needed
if args.resolutions is None:
    reslist = ["10000000", "5000000", "2500000", "1000000", "500000", "250000", "100000", "50000", "25000", "10000", "5000", "2000", "1000"]
resst = ",".join(reslist)
#Extract minimum resolution
minres = reslist[-1]

#Check for RCMC mode
if captureregion is not None:
    RCMCmode = 1

#Process ID (used to make unique sorttemp, so these are not overlapping for multiple processes in the same outdir)
uniqueid = str(uuid.uuid4())

# commands as strings
line1 = "mkdir {0}{10}sorttemp -p"
#Line 2 options
#Bowtie options
line2bowtiemask = "bowtie2 -x /mnt/md0/DataRepository/genomes/{1}/{2} --threads {3} -1 {4} -2 {5} --reorder --local --very-sensitive-local {12}{11}| pairtools parse --add-columns mapq --walks-policy mask -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"
line2bowtieall = "bowtie2 -x /mnt/md0/DataRepository/genomes/{1}/{2} --threads {3} -1 {4} -2 {5} --reorder --local --very-sensitive-local {12}{11}| pairtools parse --add-columns mapq --walks-policy all -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"
line2bowtiep2 = "bowtie2 -x /mnt/md0/DataRepository/genomes/{1}/{2} --threads {3} -1 {4} -2 {5} --reorder --local --very-sensitive-local {12}{11}| pairtools parse2 --add-columns mapq -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"
line2bowtiep2e = "bowtie2 -x /mnt/md0/DataRepository/genomes/{1}/{2} --threads {3} -1 {4} -2 {5} --reorder --local --very-sensitive-local {12}{11}| pairtools parse2 --add-columns mapq -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --expand --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"
#bwa options
line2bwamask = "bwa-mem2 mem -t {3} -SP /mnt/md0/DataRepository/genomes/{1}/{2}.fa {4} {5} {11}| pairtools parse --add-columns mapq --walks-policy mask -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"
line2bwaall = "bwa-mem2 mem -t {3} -SP /mnt/md0/DataRepository/genomes/{1}/{2}.fa {4} {5} {11}| pairtools parse --add-columns mapq --walks-policy mask -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"
line2bwap2 = "bwa-mem2 mem -t {3} -SP /mnt/md0/DataRepository/genomes/{1}/{2}.fa {4} {5} {11}| pairtools parse2 --add-columns mapq -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"
line2bwap2e = "bwa-mem2 mem -t {3} -SP /mnt/md0/DataRepository/genomes/{1}/{2}.fa {4} {5} {11}| pairtools parse2 --add-columns mapq -c /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes --assembly {2} --min-mapq {13} --drop-sam --drop-readid --expand --nproc-in {3} | pairtools sort --tmpdir {0}{10}sorttemp --nproc {3} -o {0}{6}.pairs | cat"

if aligner == "bowtie2":
    if walks == "mask":
        line2 = line2bowtiemask
    elif walks == "all":
        line2 = line2bowtieall
    elif walks == "parse2":
        line2 = line2bowtiep2
    elif walks == "expand":
        line2 = line2bowtiep2
    else:
        print("Please provide a valid option for --walks - see help.")
        exit()
elif aligner == "bwa":
    if walks == "mask":
        line2 = line2bwamask
    elif walks == "all":
        line2 = line2bwaall
    elif walks == "parse2":
        line2 = line2bwap2
    elif walks == "expand":
        line2 = line2bwap2e
    else:
        print("Please provide a valid option for --walks - see help.")
        exit()


line3 = "pairtools merge --tmpdir {0}{10}sorttemp --nproc {3} {7} | pairtools dedup --max-mismatch 1 --mark-dups --output {0}{6}.nodups.pairs.gz --output-unmapped {0}{6}.unmapped.pairs.gz --output-dups {0}{6}.dups.pairs.gz --output-stats {0}{6}.dedup.stats | cat"
line4 = "pairtools dedup --max-mismatch 1 --mark-dups --output {0}{6}.nodups.pairs.gz --output-unmapped {0}{6}.unmapped.pairs.gz --output-dups {0}{6}.dups.pairs.gz --output-stats {0}{6}.dedup.stats {7}"
line5 = "pairix {0}{6}.nodups.pairs.gz"
line6 = "bgzip -cd -@ 3 {0}{6}.nodups.pairs.gz | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly {2} /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes:{8} - {0}{6}.{8}.cool"
line7 = "cooler zoomify --nproc {3} --balance --out {0}{6}.{8}.mcool --resolutions {9} {0}{6}.{8}.cool"
line8 = "rmdir {0}{10}sorttemp"

#For running only alignment and making bams
bline1 = "mkdir {0}{6}temp -p"
bline2 = "bowtie2 -x /mnt/md0/DataRepository/genomes/{1}/{2} --threads {3} -1 {4} -2 {5} --reorder --local --very-sensitive-local | samtools view -bS -o {0}{6}.bam"
bline3 = "sambamba sort -t {3} -m 6GB --tmpdir {0}{7}temp {0}{6}.bam {0}{6}.sorted.bam && rm {0}{6}.bam"
bline4 = "sambamba merge -t {3} {0}{6}.sorted.merged.bam {7} && rm {7} {8}"
bline5 = "sambamba markdup -t {3} --tmpdir {0}{6}temp --overflow-list-size 10000000 -r {0}{6}.sorted.merged.bam {0}{6}.nodups.sorted.merged.bam && rm {0}{6}.sorted.merged.ba*"
bline6 = "sambamba markdup -t {3} --tmpdir {0}{6}temp -r {0}{6}.sorted.bam {0}{6}.nodups.sorted.bam && rm {0}{6}.sorted.ba*"
bline7 = "rmdir {0}{6}temp"

multilines = [line1, line2, line3, line5, line6, line7, line8]
lines = [line1, line2, line4, line5, line6, line7, line8]

multiblines = [bline1, bline2, bline3, bline4, bline5, bline7]
blines = [bline1, bline2, bline3, bline6, bline7]

truncmultiblines = [bline1, bline3, bline4, bline5, bline7]
truncblines = [bline1, bline3, bline6, bline7]



#Process the files depending on the run mode
if not multifile and not alignonly:
    #Include command if bams are wanted
    if keepbams:
        keepbamcmd = f"| tee >(samtools view -bS > {outdir}{outname}.bam) "
    else:
        keepbamcmd = ""
    for line in lines:
        # add file name and split by whitespace
        tokenized_line = line.format(outdir, gentype, genome, threads, file1, file2, outname, pairnamest, minres, resst, uniqueid, keepbamcmd, optbtfilter, filterexp)
        print(tokenized_line)
        # run
        sp.run(tokenized_line, shell=True, executable="/bin/bash")
elif multifile and not alignonly:
    for line in multilines:
        if line == line2:
            for x in range(len(file1)):
                #Include command if bams are wanted
                if keepbams:
                    keepbamcmd = "| tee >(samtools view -bS > {0}{1}.bam) ".format(outdir, outlist[x])
                else:
                    keepbamcmd = ""
                # add file name and split by whitespace
                tokenized_line = line.format(outdir, gentype, genome, threads, file1[x], file2[x], outlist[x], pairnamest, minres, resst, uniqueid, keepbamcmd, optbtfilter, filterexp)
                print(tokenized_line)
                # run
                sp.run(tokenized_line, shell=True, executable="/bin/bash")
        else:
            # add file name and split by whitespace
            tokenized_line = line.format(outdir, gentype, genome, threads, file1[1], file2[1], outname, pairnamest, minres, resst, uniqueid)
            print(tokenized_line)
            # run
            sp.run(tokenized_line, shell=True)
elif alignonly and multifile:
    for line in multiblines:
        if line == bline2 or line == bline3:
            for x in range(len(file1)):
                tokenized_line = line.format(outdir, gentype, genome, threads, file1[x], file2[x], outlist[x], outname)
                print(tokenized_line)
                sp.run(tokenized_line, shell=True)
        else:
            tokenized_line = line.format(outdir, gentype, genome, threads, file1[1], file2[1], outname, bamst, baist)
            print(tokenized_line)
            sp.run(tokenized_line, shell=True)
elif alignonly and not multifile:
    for line in blines:
        tokenized_line = line.format(outdir, gentype, genome, threads, file1, file2, outname, outname)
        print(tokenized_line)
        sp.run(tokenized_line, shell=True)

#After everything finishes, merge and process bams as required if doing full analysis
if keepbams and not alignonly:
    if not multifile:
        for line in truncblines:
            tokenized_line = line.format(outdir, gentype, genome, threads, file1, file2, outname, outname)
            print(tokenized_line)
            sp.run(tokenized_line, shell=True)
    elif multifile:
        for line in truncmultiblines:
            if line == bline3:
                for x in range(len(file1)):
                    tokenized_line = line.format(outdir, gentype, genome, threads, file1[x], file2[x], outlist[x], outname)
                    print(tokenized_line)
                    sp.run(tokenized_line, shell=True)
            else:
                tokenized_line = line.format(outdir, gentype, genome, threads, file1[1], file2[1], outname, bamst, baist)
                print(tokenized_line)
                sp.run(tokenized_line, shell=True)
