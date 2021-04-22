import sys
import getopt
import re

def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:o:",["infile=","outfile=",])
    except getopt.GetoptError:
        print('cutesv_inv_fliter.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('cutesv_inv_fliter.py -i <inputfile> -o <outputfile>')
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
            if re.search('^#',line1):
                outputfile.write(str(line1))
                outputfile.write('\n')
            else:
                cut2 = cut1[7]
                INFO = cut2.strip().split(';')
                sv_type = INFO[1]
                length = INFO[2].split('=')
                len = length[1]
                filter = INFO[0]
                cut3 = cut1[9]
                FORMAT = cut3.strip().split(':')
                GT = FORMAT[0]
                if GT =='0/1' or GT == '1/1':
                    if cut1[6] == 'PASS':
                        if filter == 'PRECISE':
                            if sv_type == 'SVTYPE=DUP':
                                if abs(int(len)) >= 20 and abs(int(len)) < 30000:
                                    outputfile.write(str(line1))
                                    outputfile.write('\n')
outputfile.close()
