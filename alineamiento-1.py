from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein, generic_rna
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


secuencia_epipedobates = "TTCTTTGTGTTGGTGATTTTCCTGGGCTCTTTCTATCTCATCAACCTTATCCTCGCTGTGGTGGCCATGGCATATGACGAGCAAAATGAAGCCACCATCCAGGAAGCCTTA"
result_handle= NCBIWWW.qblast("blastn", "nt",secuencia_epipedobates)
print(result_handle)

blast_records = NCBIXML.parse(result_handle)
dir(blast_records)
print(blast_records.__sizeof__)

print('1')
for b in blast_records:
    for alignment in b.alignments:
        for hsp in alignment.hsps:
            print('\n****Alineamiento****')
            print('secuencia:', alignment.title)
            print('longitud:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query[0:100] + '...')
            print(hsp.match[0:100] + '...')
            print(hsp.sbjct[0:100] + '...')
