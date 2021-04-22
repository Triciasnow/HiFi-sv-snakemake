import re
import sys
import getopt

def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:o:",["infile=","outfile=",])
    except getopt.GetoptError:
        print('syri_tra_fliter.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('syri_tra_fliter.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            inname = arg
        elif opt in ("-o", "--outfile"):
            outname = arg
    return inname,outname

if __name__ == "__main__":
    inputname,outputname = main(sys.argv[1:])
    inputfile = open(inputname)
    outputfile = open(outputname,'w')

    while 1:
        lines = inputfile.readlines(10000)
        if not lines:
            break
        for line1 in lines:
            line1 = line1.rstrip()
            if re.search('^#', line1):
                if re.search('^#C',line1):
                    outputfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'+'\n')
                    outputfile.write(line1+'\t'+'FORMAT'+'\t'+'Sample'+'\n')
                else:
                    outputfile.write(line1+'\n')
            else:
                cut1 = line1.strip().split('\t')
                chr = cut1[0]
                if chr == 'Chr10' or chr == 'Chr11' or chr == 'Chr12' or chr == 'Chr0':
                    chr1 = chr.strip('Chr')
                    cut1[0] = chr1
                    n = '\t'.join(cut1)
                else:
                    chr1 = chr.strip('Chr0')
                    cut1[0] = chr1
                    n = '\t'.join(cut1)

                if re.search('^NOTAL', cut1[2]):
                    continue
                elif re.search('^TRANS', cut1[2]):
                    if cut1[4] == '<TRANS>':
                        outputfile.write(n+'\t'+'GT'+'\t'+'./.'+'\n')
                elif re.search('^INV', cut1[2]):
                    if cut1[4] == '<INVTR>':
                        outputfile.write(n+'\t'+'GT'+'\t'+'./.'+'\n')
                else:
                    continue
outputfile.close()
