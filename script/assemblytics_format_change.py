import re
import sys
import getopt

def main(argv):
    inname=''
    outname=''
    try:
        opts,args=getopt.getopt(argv,"hi:a:o:",["infile=","anno=","outfile=",])
    except getopt.GetoptError:
        print('assemblytics_format_change.py -i <inputfile> -a <annofile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('assemblytics_format_change.py -i <inputfile> -a <annofile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            inname = arg
        elif opt in ("-a", "--afile"):
            aname = arg
        elif opt in ("-o", "--outfile"):
            outname = arg
    return inname,aname,outname

if __name__ == "__main__":
    inputname,annoname,outputname = main(sys.argv[1:])
    inputfile = open(inputname)
    annofile = open(annoname)
    outputfile = open(outputname,'w')

    while 1:
        lines = annofile.readlines(10000)
        if not lines:
            break
        for line1 in lines:
            outputfile.write(str(line1))
    annofile.close()

    id = 0
    while 2:
        lines = inputfile.readlines(10000)
        if not lines:
            break
        for line2 in lines:
            line2 = line2.rstrip()
            cut2 = line2.strip().split('\t')
            if re.search('^R', line1):
                continue
            else:
                if cut2[0] == '0':
                    continue
                else:
                    chr2 = cut2[0]
                    pos2 = cut2[1]
                    id2 = str(id)
                    id = int(id2) + 1
                    ref2 = '.'
                    qual2 = '.'
                    filter2 = 'PASS'

                    if cut2[6] == 'Deletion':
                        if abs(int(cut2[4])) >= 50 and abs(int(cut2[4])) < 200000:
                            que2 = '<DEL>'
                            svtype = 'DEL'
                            c2 = cut2[9].split(':')
                            c3 = c2[0]
                            info2 = 'SVTYPE='+svtype+';'+'END='+cut2[2]+';'+'SVLEN='+cut2[4]+';'+'CHR2='+c3
                            outputfile.write(chr2+'\t'+pos2+'\t'+id2+'\t'+ref2+'\t'+que2+'\t'+qual2+'\t'+filter2+'\t'+info2+'\t'+'GT'+'\t'+'./.'+'\n')
inputfile.close()
outputfile.close()
