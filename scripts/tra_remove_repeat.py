import sys
import getopt
import re

def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:o:",["infile=","para=","outfile=",])
    except getopt.GetoptError:
        print('tra_remove_repeat.py -i <inputfile>  -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('tra_remove_repeat.py -i <inputfile> -o <outputfile>')
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
    l = 0
    s = []
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
                pos = cut1[1]
                info = cut1[7].strip().split(';')
                type = info[3]
                end = info[6]
                pos_change1 = end.strip().split('=')
                pos_change = pos_change1[1]
                chr2 = info[5].split('=')
                chr_new = chr2[1]
                all = chr+':'+pos+'_'+chr_new+':'+pos_change
                all_back = chr_new + ':' + pos_change +'_'+ chr +':'+pos
                if all in s or all_back in s:
                    continue
                else:
                    s.append(all)
                    s.append(all_back)
                    outputfile.write(line1+'\n')
outputfile.close()
