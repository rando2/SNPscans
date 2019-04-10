import csv
window=100000
header = ""
ind={"Tame":[],"Aggr":[]}
with open("../fox5k.quantfilt.recode.vcf",'rb') as vcffile, open("./karlsonfst.quantfilt.tsv","wb") as outfile:
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
            writerbot.writerow(["CHROM","POS","KARLSSON_FST"])
        else:
            if len(line[4].split(',')) > 1: #if it's not biallelic, this formula doesn't work
                continue 
            pq={"Tame":[0.0,0.0], "Aggr":[0.0,0.0]} #p,q
            for pop in ["Tame", "Aggr"]:
                for i in ind[pop]:
                    geno = line[i].split(':')[0]
                    if geno == "0/0":
                        pq[pop][0] += 2
                    elif geno == "0/1":
                        pq[pop][0] +=1
                        pq[pop][1] +=1
                    elif geno == "1/1":
                        pq[pop][1] +=2
            if sum(pq["Tame"]) < 40 or sum(pq["Aggr"]) < 40:
                continue

            p1 = pq["Tame"][0]/sum(pq["Tame"])
            p2 = pq["Aggr"][0]/sum(pq["Aggr"])
            q1= pq["Tame"][1]/sum(pq["Tame"])
            q2= pq["Aggr"][1]/sum(pq["Aggr"])

            N = float(p1*(q2 - q1) + p2*(q1 - q2))
            D = float(p1*q2+ q1*p2)
            if D==0: #if the populations are both fixed for the same allele then Fst is 0
                Fst=0
            else:
                Fst=N/D
            writerbot.writerow(line[0:2] + [Fst])

