import os
import pandas as pd
from Bio import SeqIO
from subprocess import PIPE, run
import _RNA as RNA

input_file = 'm6a_data.txt'

fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    # cmd = 'RNAshapes -o 1 -# 1 ' + sequence
    # result = run(cmd, stdout=PIPE, universal_newlines=True, shell=True)
    # result = result.stdout
    # length = len(sequence)
    # secondary = result[length + 1: 2 * length + 1]

    secondary, mfe = RNA.fold(sequence)
    print(secondary)

    temp_file = open('temp.fa', 'w')
    temp_file.write(secondary)
    temp_file.write('\n')
    temp_file.write(sequence)
    temp_file.close()

    cmd_jar3d_IL = "java -jar jar3d/jar3d_2014-12-11.jar temp.fa jar3d/IL/3.2/lib/all.txt temp.IL.loop.txt temp.IL.sequence.txt"
    os.system(cmd_jar3d_IL)
    data_IL = pd.read_csv('temp.IL.loop.txt', index_col=1, header=0)
    data_IL = data_IL['meanCutoffScore']

    cmd_jar3d_HL = "java -jar jar3d/jar3d_2014-12-11.jar temp.fa jar3d/HL/3.2/lib/all.txt temp.HL.loop.txt temp.HL.sequence.txt"
    os.system(cmd_jar3d_HL)
    data_HL = pd.read_csv('temp.HL.loop.txt', index_col=1, header=0)
    data_HL = data_HL['meanCutoffScore']

    print(sum(data_IL>0))
    print(sum(data_HL>0))