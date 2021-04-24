import re
import sys
import getopt


def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:c:o:",["infile=","chr=","outfile="])
    except getopt.GetoptError:
        print('ref_chr.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('ref_chr.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            inname = arg
        elif opt in ("-c", "--cfile"):
            cname = arg
        elif opt in ("-o", "--outfile"):
            outname = arg
    return inname,cname,outname

if __name__ == "__main__":
    inputname,chrname,outputname = main(sys.argv[1:])
    inputfile = open(inputname)
    outputfile = open(outputname,'w')
    while 1:
        lines = inputfile.readlines(10000)
        if not lines:
            break
        for line1 in lines:
            if re.search('^>',line1):
                print(line1)
                for i in range(0,int(chrname)+1):
                    c = '>'+str(i)
                    if re.search(c,line1):
                        c1 = line1.strip('>')
                        c2 = int(c1)
                        if 0 < c2 <10:
                            c_new = '>Chr0'+str(i)
                            #print(line1)
                            outputfile.write(c_new+'\n')
                        else:
                            c_new = '>Chr'+str(i)
                            if i >= 10 or i == 0:
                                outputfile.write(c_new+'\n')
            else:
                outputfile.write(line1)
outputfile.close()
