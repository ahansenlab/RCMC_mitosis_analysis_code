#Merging aligned and parsed pairs files with or without duplicate removal, then generating cools and mcools

from sys import exit
import subprocess as sp
import argparse
import multiprocessing
import uuid

parser = argparse.ArgumentParser(description = "run pairtools and cooltools on aligned and parsed pairs files to merge them and make mcools")
parser.add_argument("--file", "-f", help = "pairs file or list of pairs files", nargs = "*")
parser.add_argument("--genome", "-g", help = "genome aligned to - mouse or human - required")
parser.add_argument("--genometype", "-y", help = "genome type - use if your genome is a modified version of a standard genome - should be one of hg19, hg38, mm10, or mm39")
parser.add_argument("--threads", "-t", help = "number of threads to use - default is 1", default = "1")
parser.add_argument("--resolutions", "-r", help = "list of resolutions to output in decreasing order - all resolutions must be a multiple of the smallest resolution - default: 10000000 5000000 2500000 1000000 500000 250000 100000 50000 25000 10000 5000 2000 1000", nargs = "*")
parser.add_argument("--out", "-o", help = "name for merged output files - defaults to name of first file")
parser.add_argument("--outdir", help = "a directory to store output files - default is current directory", default = "./")
parser.add_argument("--rmdup", "-d", help = "flag for whether or not to remove PCR duplicates after merging - set to remove if merging e.g. samples from different sequencing runs on the same library - default: not set", action = "store_true")
args = parser.parse_args()

file = args.file
genome = args.genome
gentype = args.genometype
threads = args.threads
outname = args.out
outdir = args.outdir
reslist = args.resolutions

#Check requirements are fulfilled:
condapacks = sp.run("conda list".split(), capture_output=True)
condapacksstr = str(condapacks.stdout)

if "pairtools" not in condapacksstr or "cooler" not in condapacksstr or "pairix" not in condapacksstr:
    print("Please make sure pairtools, pairix and cooler are installed in your current conda environment (check conda list)")
    exit()

#Check that files are specified
if file is None:
    print("Input files not specified - check help for formatting")
    parser.print_usage()
    exit()
elif len(file) > 1:
    filest = " ".join(file)
    if args.out is None:
        outname = file[0] + ".merged"
    else:
        outname = args.out + ".merged"
    #pairnamelist = [outdir + oname + ".pairs" for oname in outlist]
    #pairnamest = " ".join(pairnamelist)
elif len(file1) == 1:
    print("How are you going to merge one file?")
    exit()
else:
    print("Mismatch in number of input files, check arguments")
    exit()

#Check that outdir ends with a /, add one if it doesn't
if not outdir.endswith("/"):
    outdir = outdir + "/"

#Check that a sensible number of threads has been requested - more protections here are possible - at the moment users are trusted to be sensible
cpucount = multiprocessing.cpu_count()
if args.threads is None:
    print("Defaulting to one thread")
    threads = 1
elif int(args.threads) >= cpucount:
    print("Too many threads requested, resetting to default")
    threads = 1

#Check that the user has entered a valid genome to align to
if args.genome is None:
    print("Please specify a genome identifier!")
    exit()

if args.genometype is None:
    gentype = args.genome


if gentype == "mm10" or gentype == "mm39" or gentype == "mm9":
    toprint = "Aligning to mouse genome {}".format(genome)
    print(toprint)
elif gentype == "hg19" or gentype == "hg38":
    toprint = "Aligning to human genome {}".format(genome)
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

#Process ID (used to make unique sorttemp, so these are not overlapping for multiple processes in the same outdir)
uniqueid = str(uuid.uuid4())

# commands as strings
line1 = "mkdir {0}{8}sorttemp -p"
line2 = "pairtools merge --tmpdir {0}{8}sorttemp --nproc {3} -o {0}{5}.nodups.pairs.gz {4}"
line3 = "pairtools merge --tmpdir {0}{8}sorttemp --nproc {3} {4} | pairtools dedup --max-mismatch 1 --mark-dups --output {0}{5}.nodups.pairs.gz --output-unmapped {0}{5}.unmapped.pairs.gz --output-dups {0}{5}.dups.pairs.gz --output-stats {0}{5}.dedup.stats | cat"
line4 = "pairix {0}{5}.nodups.pairs.gz"
line5 = "bgzip -cd -@ 3 {0}{5}.nodups.pairs.gz | cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly {2} /mnt/md0/DataRepository/chromsizes/{1}/{2}.sorted.chrom.sizes:{6} - {0}{5}.{6}.cool"
line6 = "cooler zoomify --nproc {3} --out {0}{5}.{6}.mcool --resolutions {7} --balance {0}{5}.{6}.cool"
line7 = "rmdir {0}{8}sorttemp"

nodedup = [line1, line2, line4, line5, line6, line7]
dedup = [line1, line3, line4, line5, line6, line7]

if args.rmdup:
    print("Merging and removing PCR duplicates...")
    for line in dedup:
        # add file name and split by whitespace
        tokenized_line = line.format(outdir, gentype, genome, threads, filest, outname, minres, resst, uniqueid)
        print(tokenized_line)
        # run
        sp.run(tokenized_line, shell=True)
    print(outname, "finished!")
else:
    print("Merging...")
    for line in nodedup:
        # add file name and split by whitespace
        tokenized_line = line.format(outdir, gentype, genome, threads, filest, outname, minres, resst, uniqueid)
        print(tokenized_line)
        # run
        sp.run(tokenized_line, shell=True)
    print(outname, "finished!")
