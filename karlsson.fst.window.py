import csv
import numpy as np

window=200000
design="indiv"

header = ""
ind={"Tame":[],"Aggr":[]}
lastscaff=""
w=[]
N1=[]
N2=[]
D1=[]
D2=[]
FST1=[]
FST2=[]

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
    #eFST indicates empirical FSTs (the ones directly calculated by N/D at each site), whereas tFST is the theoretical Fst, smoothed over all sites
    num_points = float(len(Ns))
    if num_points == 0:
        print "error: no points in vectors"
        print lastscaff, w
        return 0,0
    N = np.array(Ns)
    D = np.array(Ds)
    tFST = np.sum(N) / np.sum(D) #based on Karlsson et al., Nature Genetics 2007
    meanFST= np.mean(eFST)
    s2FST = (1/num_points)*np.sum((eFST-meanFST)**2) #based on Oleksyk et al., 2008 PLoSOne
    return round(tFST,5), round(s2FST,5)

def reset_running(): #use to set a clean slate of N and D
    return [], [], []

def parse_geno(geno, tot, design="indiv"):
    #provide geno = genotype to parse (0/0 for indiv or 3,3 for pooled)
    #tot = dictionary in format [0,0] representing alleles a1 and a2
    #design = whether it should parse GT calls or DP (for deep vs shallow/pooled design)
    #assumes biallelic & GATK format (/ sep GT and , sep DP)
    if design== "indiv":
        #if geno == "0/0":
        #    tot[0] += 2
        #elif geno == "0/1":
        #    tot[0] +=1
        #    tot[1] +=1
        #elif geno == "1/1":
        #    tot[1] +=2
        
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

with open("/home/rando2/deepTA/vcf/fox5k.quantfilt.recode.vcf",'rb') as vcffile, open("./karlsonfst.quantfilt.window.tsv","wb") as outfile:
    vcf = csv.reader(vcffile, delimiter='\t')
    writerbot = csv.writer(outfile, delimiter='\t')
    for line in vcf:
        if line[0][0:2] == "##":
            continue
        if header == "":
            header = line #['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',...
            for i in range(9,len(header)):
                pop=header[i].split('_')[1]
                if pop != "F1":
                    ind[pop].append(i)
            writerbot.writerow(["CHROM","STARTPOS","ENDPOS","NUMSNP","KARLSSON_FST","FST_VAR"])
        else:
            scaff = line[0]
            pos = int(line[1])
            if len(line[4].split(',')) > 1: #if it's not biallelic, this formula doesn't work
                continue 
            if scaff != lastscaff:
                if lastscaff != "":
                    if w[0] + window < scafflen[scaff]: #this shouldn't really happen, but just in case
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

            if pos > w[0] + window: #if the position is outside the window and there is a next window
                FST, s2FST = calc_metrics(N1, D1, FST1) #calculate the metrics
                writerbot.writerow([scaff,w[0],w[0]+window, len(N1),FST, s2FST]) #record them
                w[0] += window #increment onto next window
                N1, D1, FST1 = reset_running()
            if pos > w[1] + window:
                FST, s2FST = calc_metrics(N2, D2, FST2)
                writerbot.writerow([scaff,w[1],w[1]+window, len(N2),FST, s2FST])
                w[1] += window
                N2, D2,FST2 = reset_running()

            pq={"Tame":[0.0,0.0], "Aggr":[0.0,0.0]} #p,q
            for pop in ["Tame", "Aggr"]: #tabulate the alleles for each of them
                for i in ind[pop]: #the values in the dict are NOT the names, they are the indices in the header 
                    geno = line[i].split(':')[0]
                    pq[pop]=parse_geno(geno, pq[pop], "indiv")

            if sum(pq["Tame"]) < 40 or sum(pq["Aggr"]) < 40: #if fewer than 20 animals genotyped per pop, skip
                continue
            n1=sum(pq["Tame"]) #n is total number of alleles
            n2=sum(pq["Aggr"]) #n is total number of alleles

            p1 = pq["Tame"][0]/n1
            p2 = pq["Aggr"][0]/n2
            q1= pq["Tame"][1]/n1
            q2= pq["Aggr"][1]/n2
            
            #This formula is definitely not correct (assumes = sample sizes), but it's the one Jen told me to use to match her analysis
            #N = float(p1*(q2 - q1) + p2*(q1 - q2))
            #D = float(p1*q2+ q1*p2)

            h1= (pq["Tame"][0]*pq["Tame"][1]) / (n1*(n1-1))
            h2= (pq["Aggr"][0]*pq["Aggr"][1]) / (n2*(n2-1))

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
                print "hopefully", pos, "is near the end of",scaff, "at", scafflen[scaff]
            if pos >= w[1] & pos < w[1] + window:
                N2.append(N)
                D2.append(D)
                FST2.append(FST)
            else:
                print "problem with windows", pos, scaff
    FST, s2FST = calc_metrics(N2, D2, FST2)
    writerbot.writerow([lastscaff,w[1],w[1]+window, len(N2),FST, s2FST])
