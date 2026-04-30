import re
import gzip
import logging
from time import localtime, strftime

timeformat = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(level=logging.INFO, format=timeformat)

SEED1_LE = "\w?AATCTGTCACCGACGACAGATAATTTGTCACTGT|T\w?ATCTGTCACCGACGACAGATAATTTGTCACTGT|TA\w?TCTGTCACCGACGACAGATAATTTGTCACTGT|TAA\w?CTGTCACCGACGACAGATAATTTGTCACTGT|TAAT\w?TGTCACCGACGACAGATAATTTGTCACTGT|TAATC\w?GTCACCGACGACAGATAATTTGTCACTGT|TAATCT\w?TCACCGACGACAGATAATTTGTCACTGT|TAATCTG\w?CACCGACGACAGATAATTTGTCACTGT|TAATCTGT\w?ACCGACGACAGATAATTTGTCACTGT|TAATCTGTC\w?CCGACGACAGATAATTTGTCACTGT|TAATCTGTCA\w?CGACGACAGATAATTTGTCACTGT|TAATCTGTCAC\w?GACGACAGATAATTTGTCACTGT|TAATCTGTCACC\w?ACGACAGATAATTTGTCACTGT|TAATCTGTCACCG\w?CGACAGATAATTTGTCACTGT|TAATCTGTCACCGA\w?GACAGATAATTTGTCACTGT|TAATCTGTCACCGAC\w?ACAGATAATTTGTCACTGT|TAATCTGTCACCGACG\w?CAGATAATTTGTCACTGT|TAATCTGTCACCGACGA\w?AGATAATTTGTCACTGT|TAATCTGTCACCGACGAC\w?GATAATTTGTCACTGT|TAATCTGTCACCGACGACA\w?ATAATTTGTCACTGT|TAATCTGTCACCGACGACAG\w?TAATTTGTCACTGT|TAATCTGTCACCGACGACAGA\w?AATTTGTCACTGT|TAATCTGTCACCGACGACAGAT\w?ATTTGTCACTGT|TAATCTGTCACCGACGACAGATA\w?TTTGTCACTGT|TAATCTGTCACCGACGACAGATAA\w?TTGTCACTGT|TAATCTGTCACCGACGACAGATAAT\w?TGTCACTGT|TAATCTGTCACCGACGACAGATAATT\w?GTCACTGT|TAATCTGTCACCGACGACAGATAATTT\w?TCACTGT|TAATCTGTCACCGACGACAGATAATTTG\w?CACTGT|TAATCTGTCACCGACGACAGATAATTTGT\w?ACTGT|TAATCTGTCACCGACGACAGATAATTTGTC\w?CTGT|TAATCTGTCACCGACGACAGATAATTTGTCA\w?TGT|TAATCTGTCACCGACGACAGATAATTTGTCAC\w?GT|TAATCTGTCACCGACGACAGATAATTTGTCACT\w?T|TAATCTGTCACCGACGACAGATAATTTGTCACTG\w?"
SEED1_RE = "\w?AATAATTTGTCACAACGACATATAATTAGTCACT|C\w?ATAATTTGTCACAACGACATATAATTAGTCACT|CA\w?TAATTTGTCACAACGACATATAATTAGTCACT|CAA\w?AATTTGTCACAACGACATATAATTAGTCACT|CAAT\w?ATTTGTCACAACGACATATAATTAGTCACT|CAATA\w?TTTGTCACAACGACATATAATTAGTCACT|CAATAA\w?TTGTCACAACGACATATAATTAGTCACT|CAATAAT\w?TGTCACAACGACATATAATTAGTCACT|CAATAATT\w?GTCACAACGACATATAATTAGTCACT|CAATAATTT\w?TCACAACGACATATAATTAGTCACT|CAATAATTTG\w?CACAACGACATATAATTAGTCACT|CAATAATTTGT\w?ACAACGACATATAATTAGTCACT|CAATAATTTGTC\w?CAACGACATATAATTAGTCACT|CAATAATTTGTCA\w?AACGACATATAATTAGTCACT|CAATAATTTGTCAC\w?ACGACATATAATTAGTCACT|CAATAATTTGTCACA\w?CGACATATAATTAGTCACT|CAATAATTTGTCACAA\w?GACATATAATTAGTCACT|CAATAATTTGTCACAAC\w?ACATATAATTAGTCACT|CAATAATTTGTCACAACG\w?CATATAATTAGTCACT|CAATAATTTGTCACAACGA\w?ATATAATTAGTCACT|CAATAATTTGTCACAACGAC\w?TATAATTAGTCACT|CAATAATTTGTCACAACGACA\w?ATAATTAGTCACT|CAATAATTTGTCACAACGACAT\w?TAATTAGTCACT|CAATAATTTGTCACAACGACATA\w?AATTAGTCACT|CAATAATTTGTCACAACGACATAT\w?ATTAGTCACT|CAATAATTTGTCACAACGACATATA\w?TTAGTCACT|CAATAATTTGTCACAACGACATATAA\w?TAGTCACT|CAATAATTTGTCACAACGACATATAAT\w?AGTCACT|CAATAATTTGTCACAACGACATATAATT\w?GTCACT|CAATAATTTGTCACAACGACATATAATTA\w?TCACT|CAATAATTTGTCACAACGACATATAATTAG\w?CACT|CAATAATTTGTCACAACGACATATAATTAGT\w?ACT|CAATAATTTGTCACAACGACATATAATTAGTC\w?CT|CAATAATTTGTCACAACGACATATAATTAGTCA\w?T|CAATAATTTGTCACAACGACATATAATTAGTCAC\w?"
SEED2 = "AA\w?ATCAGCTTCCAGCTGCCTT|AA[GCT]TATCAGCTTCCAGCTGCCTT|AAT\w?TCAGCTTCCAGCTGCCTT|AAT\w?ATCAGCTTCCAGCTGCCTT|AATA\w?CAGCTTCCAGCTGCCTT|AATA\w?TCAGCTTCCAGCTGCCTT|AATAT\w?AGCTTCCAGCTGCCTT|AATAT\w?CAGCTTCCAGCTGCCTT|AATATC\w?GCTTCCAGCTGCCTT|AATATC\w?AGCTTCCAGCTGCCTT|AATATCA\w?CTTCCAGCTGCCTT|AATATCA\w?GCTTCCAGCTGCCTT|AATATCAG\w?TTCCAGCTGCCTT|AATATCAG\w?CTTCCAGCTGCCTT|AATATCAGC\w?TCCAGCTGCCTT|AATATCAGC\w?TTCCAGCTGCCTT|AATATCAGCT\w?CCAGCTGCCTT|AATATCAGCT\w?TCCAGCTGCCTT|AATATCAGCTT\w?CAGCTGCCTT|AATATCAGCTT\w?CCAGCTGCCTT|AATATCAGCTTC\w?AGCTGCCTT|AATATCAGCTTC\w?CAGCTGCCTT|AATATCAGCTTCC\w?GCTGCCTT|AATATCAGCTTCC\w?AGCTGCCTT|AATATCAGCTTCCA\w?CTGCCTT|AATATCAGCTTCCA\w?GCTGCCTT|AATATCAGCTTCCAG\w?TGCCTT|AATATCAGCTTCCAG\w?CTGCCTT|AATATCAGCTTCCAGC\w?GCCTT|AATATCAGCTTCCAGC\w?TGCCTT|AATATCAGCTTCCAGCT\w?CCTT|AATATCAGCTTCCAGCT\w?GCCTT|AATATCAGCTTCCAGCTG\w?CTT|AATATCAGCTTCCAGCTG\w?CCTT|AATATCAGCTTCCAGCTGC\w?TT|AATATCAGCTTCCAGCTGC\w?CTT|AATATCAGCTTCCAGCTGCC\w?T|AATATCAGCTTCCAGCTGCC[AGC]TT|AATATCAGCTTCCAGCTGCCT\w?"


def filter_fq(fq1, fq2, biotin_primer_base, umi_length, seed1, seed2):
    """
    过滤FASTQ文件，基于seed1和seed2模式匹配
    
    参数:
    fq1: 第一个FASTQ文件路径
    fq2: 第二个FASTQ文件路径
    biotin_primer_base: 生物素引物使用的cas12k_handle碱基长度
    umi_length: UMI的长度
    """

    startTime = strftime(timeformat, localtime())
    logging.info(f"@@ filter fastq: {fq1}, {fq2}")

    out_pass_fq1 = f'FliterPass_{fq1}'
    out_pass_fq2 = f'FliterPass_{fq2}'
    out_fail_fq1 = f'FliterFail_{fq1}'
    out_fail_fq2 = f'FliterFail_{fq2}'

    line_number = 0
    pass_line_number = 0
    fail_line_number = 0

    reads1_list, reads2_list = list(), list()

    with (
        gzip.open(fq1, 'rt') as f1,
        gzip.open(fq2, 'rt') as f2,
        gzip.open(out_pass_fq1, 'wt') as pass_out1,
        gzip.open(out_pass_fq2, 'wt') as pass_out2,
        gzip.open(out_fail_fq1, 'wt') as fail_out1,
        gzip.open(out_fail_fq2, 'wt') as fail_out2
    ):
        
        for line1, line2 in zip(f1, f2):
            line_number += 1
            reads1_list.append(line1.strip())
            reads2_list.append(line2.strip())

            if line_number%4 == 2 and len(reads1_list) > 4:
                tmp_r1 = reads1_list[1]
                tmp_r2 = reads2_list[1]

                if re.search(seed1 ,tmp_r1[:biotin_primer_base + 15]) and re.search(seed2 ,tmp_r2[:35]):
                    # R1-R2
                    pass
                    umi = tmp_r2[:umi_length]
                    name = f'{reads1_list[0].split(" ")[0]}_{umi}'
                    read1_name = f'{name} {" ".join(reads1_list[0].split(" ")[1:])}'
                    read2_name = f'{name} {" ".join(reads2_list[0].split(" ")[1:])}'
                    reads1_list[0] = read1_name
                    reads2_list[0] = read2_name

                    pass_out1.write("\n".join(reads1_list[:4]) + "\n")
                    pass_out2.write("\n".join(reads2_list[:4]) + "\n")
                    pass_line_number += 1
                    
                elif re.search(seed1 ,tmp_r2[:biotin_primer_base + 15]) and re.search(seed2 ,tmp_r1[:35]):
                    # R2-R1
                    pass
                    umi = tmp_r1[:umi_length]
                    name = f'{reads1_list[0].split(" ")[0]}_{umi}'
                    # 交换R1和R2名字
                    read1_name = f'{name} {" ".join(reads2_list[0].split(" ")[1:])}'
                    read2_name = f'{name} {" ".join(reads1_list[0].split(" ")[1:])}'
                    reads1_list[0] = read1_name
                    reads2_list[0] = read2_name

                    pass_out1.write("\n".join(reads2_list[:4]) + "\n")
                    pass_out2.write("\n".join(reads1_list[:4]) + "\n")
                    pass_line_number += 1

                else:
                    # 没有找到cas12k handle或者UMI
                    pass
                    fail_out1.write("\n".join(reads1_list[:4]) + "\n")
                    fail_out2.write("\n".join(reads2_list[:4]) + "\n")
                    fail_line_number += 1

                reads1_list = reads1_list[:4]
                reads2_list = reads2_list[:4]

        # 最后一个read
        tmp_r1 = reads1_list[1]
        tmp_r2 = reads2_list[1]

        if re.search(seed1 ,tmp_r1[:biotin_primer_base + 15]) and re.search(seed2 ,tmp_r2[:35]):
            # R1-R2
            pass
            umi = tmp_r2[:umi_length]
            name = f'{reads1_list[0].split(" ")[0]}_{umi}'
            read1_name = f'{name} {" ".join(reads1_list[0].split(" ")[1:])}'
            read2_name = f'{name} {" ".join(reads2_list[0].split(" ")[1:])}'
            reads1_list[0] = read1_name
            reads2_list[0] = read2_name

            pass_out1.write("\n".join(reads1_list[:4]) + "\n")
            pass_out2.write("\n".join(reads2_list[:4]) + "\n")
            pass_line_number += 1
            
        elif re.search(seed1 ,tmp_r2[:biotin_primer_base + 15]) and re.search(seed2 ,tmp_r1[:35]):
            # R2-R1
            pass
            umi = tmp_r1[:umi_length]
            name = f'{reads1_list[0].split(" ")[0]}_{umi}'
            # 交换R1和R2名字
            read1_name = f'{name} {" ".join(reads2_list[0].split(" ")[1:])}'
            read2_name = f'{name} {" ".join(reads1_list[0].split(" ")[1:])}'
            reads1_list[0] = read1_name
            reads2_list[0] = read2_name

            pass_out1.write("\n".join(reads2_list[:4]) + "\n")
            pass_out2.write("\n".join(reads1_list[:4]) + "\n")
            pass_line_number += 1

        else:
            # 没有找到cas12k handle或者UMI
            pass
            fail_out1.write("\n".join(reads1_list[:4]) + "\n")
            fail_out2.write("\n".join(reads2_list[:4]) + "\n")
            fail_line_number += 1

        
    endTime = strftime(timeformat, localtime())
    logging.info(f"@@ Filter Pass {pass_line_number} fastq; Fail {fail_line_number} fastq")
    logging.info(f"@@ filter fastq end -- Run time : {endTime} - {startTime}")
    
    return out_pass_fq1, out_pass_fq2, out_fail_fq1, out_fail_fq2


if __name__ == "__main__":
    """
    主函数，处理命令行参数并调用filter_fq函数
    """
    parser = argparse.ArgumentParser(description='Filter FASTQ files based on seed patterns')
    
    # 必需参数
    parser.add_argument('--fq1', help='First FASTQ file path')
    parser.add_argument('--fq2', help='Second FASTQ file path')
    parser.add_argument('--mode', help='Library mode, LE or RE')
    
    # 可选参数
    parser.add_argument('--biotin_primer_base', type=int, default=70, 
                      help='Length of cas12k_handle for biotin primer (default: 70)')
    parser.add_argument('--umi_length', type=int, default=8, 
                      help='Length of UMI (default: 8)')
    
    args = parser.parse_args()

    if args.mode == 'LE':
        seed1 = SEED1_LE
        seed2 = SEED2
    elif args.mode == 'RE':
        seed1 = SEED1_RE
        seed2 = SEED2
    else:
        parser.error("Unknown library mode. Please use LE or RE.")

    
    # 调用过滤函数
    pass_fq1, pass_fq2, fail_fq1, fail_fq2 = filter_fq(
        args.fq1, 
        args.fq2, 
        args.biotin_primer_base, 
        args.umi_length,
        seed1,
        seed2
    )
    
    logging.info(f"@@ Output files:")
    logging.info(f"@@ Passed reads: {pass_fq1}, {pass_fq2}")
    logging.info(f"@@ Failed reads: {fail_fq1}, {fail_fq2}")