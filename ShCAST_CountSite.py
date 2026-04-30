import sys,os
import logging
from time import localtime, strftime
from collections import Counter
from optparse import OptionParser as OP

timeformat = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(level=logging.INFO, format=timeformat)


def IS_analysis(samfilename):
    if os.path.exists(samfilename) :
        pass
    else:
        logging.error("Error: %s Files not exists !"%samfilename)
        sys.exit(0)

    #--- IS_analysis ---
    # startTime = strftime(timeformat, localtime())
    logging.info("@@ IS_analysis %s",samfilename)
    reads1_number , proper_number = 0,0

    out_filename = f"insertion_site_{samfilename.replace('.sam','')}.txt"
    
    with open(samfilename) as f, open(out_filename,"w") as fh:
        for line in f:
            if line.startswith("@"):
                continue
            lineL = line.split("\t")
            reads1_number += 1
            r1_strat = int(lineL[3])
            r2_strat = int(lineL[7])
            if abs(r2_strat - r1_strat) > 10000:
                continue

            # the mate is mapped to the reverse strand #flag include 32
            if bin(int(lineL[1]))[-6] == "1":
                if (r2_strat - r1_strat) < -10:
                    continue
                pos = int(lineL[3])
                fh.write(f"{lineL[2]}\t{pos}\t+\t")


            # the read is mapped to the reverse strand #flag include 16
            elif bin(int(lineL[1]))[-5] == "1":
                if (r2_strat - r1_strat) > 10:
                    continue
                pos = int(lineL[3]) + len(lineL[9])
                fh.write(f"{lineL[2]}\t{pos}\t+\t")
            else:
                continue
                #print bin(int(lineL[1])),line,
            umi = lineL[0].split("_")[-1]

            fh.write(f"{umi}\n")
            proper_number += 1
    
    os.system(f'grep "-" {out_filename} | sort -k 2 -n -o tmp.neg.txt ')
    os.system(f'grep "+" {out_filename} | sort -k 2 -n -o tmp.pos.txt ')
    os.system(f'cat tmp.neg.txt tmp.pos.txt > sort.{out_filename}')
    os.system('rm -f tmp.neg.txt tmp.pos.txt')
    logging.info(f"------mapping file has {reads1_number} reads1, {proper_number} is proper reads------")
    # endTime = strftime(timeformat, localtime())
    return f"sort.{out_filename}"


def count(resultfile, outfile):
    if os.path.exists(resultfile) :
        pass
    else:
        logging.error("Error: %s Files not exists !"%resultfile)
        sys.exit(0)

    #--- IS_count ---
    logging.info("@@ IS_count %s",resultfile)
    
    countL = list()
    countS = set()
    prepos = 0
    prestrand = ""

    with open(resultfile) as f, open(outfile,"w") as fh:
        fh.write("#Chr\tPos\tstrand\tumi_count\treads_count\n")

        for line in f:
            lineL = line.strip().split('\t')
            pos = int(lineL[1])
            strand = lineL[2]
            nameL = lineL[:-1]
            umi = lineL[-1]

            if strand != prestrand and countL:
                reads_count = len(countL)
                umi_count = len(countS)
                line_result = "%s\t%d\t%d\n" %(Counter(countL).most_common(1)[0][0],umi_count,reads_count)
                fh.write(line_result)
                countL = list()
                countS = set()

            if pos - prepos > 6 and countL:
                reads_count = len(countL)
                umi_count = len(countS)
                line_result = "%s\t%d\t%d\n" %(Counter(countL).most_common(1)[0][0],umi_count,reads_count)
                fh.write(line_result)
                countL = list()
                countS = set()

            prepos = pos
            prestrand = strand
            countL.append("\t".join(nameL))
            countS.add(umi)

        # 最后一个site
        line_result = "%s\t%d\t%d\n" %(Counter(countL).most_common(1)[0][0],umi_count,reads_count)
        fh.write(line_result)


    return

if __name__ == "__main__":
    """
    主函数，处理命令行参数
    """
    parser = argparse.ArgumentParser(description='Count insertion sites')
    parser.add_argument('--samfilename', type=str, required=True, help='Input SAM file')
    parser.add_argument('--outfile', type=str, required=True, help='Output file')
    args = parser.parse_args()

    startTime = strftime(timeformat, localtime())
    
    IS_analysis(args.samfilename)
    count(IS_analysis(args.samfilename), args.outfile)
    
    endTime = strftime(timeformat, localtime())
    logging.info(f"count result_file -- Run time : %s - %s"%(startTime, endTime))
    