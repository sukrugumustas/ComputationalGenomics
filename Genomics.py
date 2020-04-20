#!usr/bin/python
import re
import csv

tab = str.maketrans("ACTG", "TGAC")


def main():
    sequences = {"NEPAL": get_and_validate_input("nepal_COVID19.txt"),
                 "WUHAN": get_and_validate_input("wuhan_COVID19.txt")}
    mer_array = {}
    for origin in sequences:
        sequence = sequences[origin]
        len_seq = len(sequence)
        print('Origin [%s] -> Output file: %s_output.csv' % (origin, origin))
        file = open('%s_output.csv' % origin, 'w')
        csv_writer = csv.writer(file)
        csv_writer.writerow(['Length', 'K-Mer', 'Count', 'Revcomp', 'Count'])
        for k in range(5, len_seq):
            for i in range(len_seq - k + 1):
                if sequence[i: i + k] in mer_array:
                    mer_array[sequence[i:i + k]] += 1
                else:
                    mer_array[sequence[i: i + k]] = 1
            for i in list(mer_array.keys()):
                if mer_array[i] <= 1:
                    del mer_array[i]
            mer_array = {k: v for k, v in sorted(mer_array.items(), key=lambda item: item[1], reverse=True)}
            if len(mer_array) > 0:
                for kmer in mer_array:
                    revcomp = rev_comp(kmer)
                    rev_comp_val = 0
                    if mer_array.get(revcomp) is not None:
                        rev_comp_val = mer_array[revcomp]
                    csv_writer.writerow([k, kmer, mer_array[kmer], revcomp, rev_comp_val])
            mer_array.clear()
        file.close()


def get_and_validate_input(filename):
    file = open("./res/" + filename, "r")
    lines = file.read()
    file.close()
    lines = lines.upper()
    lines = re.sub(r'[\s]', "", lines)
    if not re.match(r'[AGCT]{5,}', lines):
        raise ValueError("[" + filename + "] Input is invalid! [Must only contain A, G, C, T and be longer than 5!]")
    return lines


def rev_comp(sub_sequence):
    return sub_sequence.translate(tab)[::-1]


if __name__ == '__main__':
    main()
