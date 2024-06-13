import random
import numpy


UMI_LEN = 9
NON_T_TAIL = 2
NUCLS = ['A', 'C', 'G', 'T']
VNUCLS =  ['A', 'C', 'G']


def get_random_seq(length, non_t_tail=0):
    seq = ""
    for i in range(length-non_t_tail):
        index = random.randint(0, len(NUCLS) - 1)
        seq += NUCLS[index]
    for i in range(non_t_tail):
        index = random.randint(0, len(VNUCLS) - 1)
        seq += VNUCLS[index]
    return seq


umi_set = set()
while len(umi_set) < 40:
    umi = get_random_seq(UMI_LEN, NON_T_TAIL)
    if umi in umi_set:
        continue
    umi_set.add(umi)
    count = 1 + round(numpy.random.gamma(0.5, 10, 1)[0])
    for i in range(count):
        print(umi)