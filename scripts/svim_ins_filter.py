import sys
import getopt
import re

def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:q:o:",["infile=","para=","outfile="])
    except getopt.GetoptError:
        print('svim_ins_filter.py -i <inputfile> -q <quality> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('svim_ins_filter.py -i <inputfile> -q <quality> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            inname = arg
        elif opt in ("-q", "--qfile"):
            qname = arg
        elif opt in ("-o", "--outfile"):
            outname = arg
    return inname,qname,outname

if __name__ == "__main__":
    inputname,qualname,outputname = main(sys.argv[1:])
    inputfile = open(inputname)
    outputfile = open(outputname,'w')

    while 1:
        lines = inputfile.readlines(10000)
        if not lines:
            break
        for line1 in lines:
            line1 = line1.rstrip()
            cut1 = line1.strip().split('\t')
            if re.search('^#',line1):
                outputfile.write(str(line1))
                outputfile.write('\n')
            else:
                chr = cut1[0]
                qual = cut1[5]
                filter = cut1[6]
                cut2 = cut1[7]
                INFO = cut2.strip().split(';')
                sv_type = INFO[0]
                length = INFO[2].split('=')
                len = length[1]
                cut3 = cut1[9]
                FORMAT = cut3.strip().split(':')
                GT = FORMAT[0]
                if chr != '0':
                    if GT == '0/1' or GT == '1/1':
                        if filter == 'PASS':
                            if int(qual) >= int(qualname):
                                if sv_type == 'SVTYPE=INS':
                                    if abs(int(len)) >= 50 and abs(int(len)) < 200000:
                                        outputfile.write(str(line1))
                                        outputfile.write('\n')

outputfile.close()
