import csv
import numpy as np

#Change these for different analyses
window=500000
design="pooled" #"pooled"
populations= ["tame", "aggr"] #["wild","farm"]

#Set up global variables
header = ""
populations = [x.title() for x in populations]
ind=dict(zip(populations,[[] for x in populations]))

if design == "pooled":
    min_alleles = 10
elif design == "indiv":
    min_alleles = 40
else:
    print "design not recognized"
    exit(1)
lastscaff=""
w=[]
N1=[]
N2=[]
D1=[]
D2=[]
FST1=[]
FST2=[]
total_sites = 0

scafflen = {}
with open("/home/rando2/wild/FULL_INDEX_SCAFFOLD_LENGTHS.csv",'rb') as lenfile:
    sclens = csv.reader(lenfile, delimiter=',')
    for line in sclens:
        scnum, scaff, sclen = line
        if int(scnum) > 3154:
            break
        else:
            scafflen[scaff] = int(sclen)

def calc_metrics(Ns,Ds, eFST):
    num_points = float(len(Ns))
    if num_points == 0:
        return 0,0
    else:
        N = np.array(Ns)
        D = np.array(Ds)
        tFST = np.sum(N) / np.sum(D) #based on Karlsson et al., Nature Genetics 2007
        meanFST= np.mean(eFST)
        s2FST = (1/num_points)*np.sum((eFST-meanFST)**2) #based on Oleksyk et al., 2008 PLoSOne
        return round(tFST,5) , round(s2FST,5)

def calc_var(FST, meanFST):
    num_points = float(len(Ns))
    s2FST = (1/num_points)*np.sum((FST-meanFST)**2) #based on Oleksyk et al., 2008 PLoSOne

def reset_running(): #use to set a clean slate of N and D
    return [], [], []

def parse_geno(geno, tot, design="indiv"):
    #provide geno = genotype to parse (0/0 for indiv or 3,3 for pooled)
    #tot = dictionary in format [0,0] representing alleles a1 and a2
    #design = whether it should parse GT calls or DP (for deep vs shallow/pooled design)
    #assumes biallelic & GATK format (/ sep GT and , sep DP)
    if '.' not in geno:
        if design== "indiv":
            geno=geno.split('/')
            for i in geno:
                tot[int(i)] +=1 #add the number of alleles
        elif design == "pooled":
                geno = geno.split(',')
                for i in [0,1]:
                    tot[i] += int(geno[i]) #add the number of reads
        else:
            print "design of study not specified correctly"
            exit(1)
    return tot

if design == "indiv":
    infile="/home/rando2/deepTA/vcf/fox5k.quantfilt.recode.vcf"
    outfile="./karlsonfst.quantfilt.window" + str(window) + ".tsv"
if design == "pooled":
    infile="/home/rando2/wild/2018/vv2align/vcf/wild.vv22.quantfilt.sort.recode.vcf"
    outfile="/home/rando2/wild/2018/vv2align/fst/karlsonfst.quantfilt.window" + str(window) + "." + populations[0][0] + populations[1][0] + ".tsv"
    print outfile

with open(infile,'rb') as vcffile, open(outfile,"wb") as outfile:
    vcf = csv.reader(vcffile, delimiter='\t')
    writerbot = csv.writer(outfile, delimiter='\t')
    for line in vcf:
        if line[0][0:2] == "##":
            continue
        if header == "":
            header = line #['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',...
            for i in range(9,len(header)):
                if '_' in header[i]: #idiosyncratic features of the "sample" names in my dataset- this is for the deep data
                    pop=header[i].split('_')[1]
                else: #the sample names aren't capitalized in the wild ata
                    pop=header[i].title()
                if pop in populations:
                    ind[pop].append(i)
                if pop in ["Tame","Aggr","Conv"] and "Farm" in populations:
                    ind["Farm"].append(i)
                if pop in ["Maryland","Newf"] and "Wild" in populations:
                    ind["Wild"].append(i)
            writerbot.writerow(["CHROM","STARTPOS","ENDPOS","NUMSNP","KARLSSON_FST","FST_VAR"])
        else:
            scaff = line[0]
            pos = int(line[1])
            if len(line[4].split(',')) > 1: #if it's not biallelic, this formula doesn't work
                continue 
            total_sites +=1
            if scaff != lastscaff:
                print scaff
                if lastscaff != "":
                    if w[0] + window < scafflen[scaff] or w[0] == 0: #first case shouldn't happen, second may affect short scaffolds
                        FST, s2FST = calc_metrics(N1, D1, FST1)
                        writerbot.writerow([lastscaff,w[0],w[0]+window, len(N1),FST, s2FST])
                    FST, s2FST = calc_metrics(N2, D2, FST2)
                    writerbot.writerow([lastscaff,w[1],w[1]+window, len(N2),FST, s2FST])
                #set start for next scaff
                w=[0,window/2] 
                lastscaff=scaff
                N1, D1, FST1 = reset_running()
                N2, D2, FST2 = reset_running()
                #continue

            while pos > w[0] + window: #if the position is outside the window
                FST, s2FST = calc_metrics(N1, D1, FST1) #calculate the metrics
                writerbot.writerow([scaff,w[0],w[0]+window, len(N1),FST, s2FST]) #record them
                w[0] += window #increment onto next window
                N1, D1, FST1 = reset_running()
            while pos > w[1] + window:
                FST, s2FST = calc_metrics(N2, D2, FST2)
                writerbot.writerow([scaff,w[1],w[1]+window, len(N2),FST, s2FST])
                w[1] += window
                N2, D2,FST2 = reset_running()

            pq=dict(zip(populations,[[0.0,0.0] for x in populations])) #p,q
            for p in populations: #tabulate the alleles for each of them
                for i in ind[p]: #the values in the dict are NOT the names, they are the indices in the header 
                    gen_dict = dict(zip(line[8].split(':'),line[i].split(':')))
                    if design == "indiv":
                        geno = gen_dict['GT'] 
                    else:
                        geno = gen_dict['AD']
                    pq[p]=parse_geno(geno, pq[p], design)
            #print line
            #print pq

            if sum(pq[populations[0]]) < min_alleles or sum(pq[populations[1]]) < min_alleles: #if fewer than 20 animals genotyped per pop, skip
                continue
            n1=sum(pq[populations[0]]) #n is total number of alleles
            n2=sum(pq[populations[1]]) #n is total number of alleles

            p1 = pq[populations[0]][0]/n1
            p2 = pq[populations[1]][0]/n2
            q1= pq[populations[0]][1]/n1
            q2= pq[populations[1]][1]/n2
            
            #This formula is definitely not correct (assumes = sample sizes), but it's the one Jen told me to use to match her analysis
            #N = float(p1*(q2 - q1) + p2*(q1 - q2))
            #D = float(p1*q2+ q1*p2)

            h1= (pq[populations[0]][0]*pq[populations[0]][1]) / (n1*(n1-1))
            h2= (pq[populations[1]][0]*pq[populations[1]][1]) / (n2*(n2-1))

            N=(p1 - p2)**2 - h1/n1- h2/n2
            D= N + h1 + h2

            if D == 0:
                FST = 0
            else:
                FST = N/D
            
            if pos >= w[0] & pos < w[0] + window:
                N1.append(N)
                D1.append(D)
                FST1.append(FST)
            else:
                print pos, "doesn't fit in w0", w[0], "on", scaff, "which is this long:", scafflen[scaff]
            if pos >= w[1] & pos < w[1] + window:
                N2.append(N)
                D2.append(D)
                FST2.append(FST)
            else:
                print pos, "doesn't fit in w1", w[1], "on", scaff, "which is this long:", scafflen[scaff]
    if w[0] + window < scafflen[scaff] or w[0] == 0:
        FST, s2FST = calc_metrics(N1, D1, FST1)
        writerbot.writerow([lastscaff,w[0],w[0]+window, len(N1),FST, s2FST])
    FST, s2FST = calc_metrics(N2, D2, FST2)
    writerbot.writerow([lastscaff,w[1],w[1]+window, len(N2),FST, s2FST])
    print design, populations, total_sites
