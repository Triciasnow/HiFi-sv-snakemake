import re
import sys
import getopt

def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:o:",["infile=","outfile=",])
    except getopt.GetoptError:
        print('vcf_fliter.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('vcf_fliter.py -i <inputfile> -o <outputfile>')
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
            cut1 = line1.strip().split('\t')
            if re.search('^#', line1):
                if re.search('^#C', line1):
                    outputfile.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO    FORMAT  VF36_Assemblytics'+'\n')
                else:
                    outputfile.write(line1+'\n')
            else:
                outputfile.write(line1+'\t'+'GT'+'\t'+'./.'+'\n')
