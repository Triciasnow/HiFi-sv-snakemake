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
            if re.search('^>',line1):
                #print(line1)
                if re.search('>ChrUn',line1):
                    outputfile.write('>Chr0'+'\n')
                else:
                    n = str(line1).strip('\n')
                    outputfile.write(n+'\n')
            else:
                outputfile.write(line1)
outputfile.close()
