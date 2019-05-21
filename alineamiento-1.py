from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein, generic_rna
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys

output_filename = "output.txt"
f = open(output_filename,'w')
sys.oudout = f
epipedobates_anthonyi = SeqIO.read("epipedobates-anthonyi-parcial.fasta", "fasta")
print("\n----Descripcion del gen Epipedobates anthonyi----")
print(epipedobates_anthonyi)
print("\nSecuencia completa del gen :")
print(epipedobates_anthonyi.seq)

epipedobates_boulengei = SeqIO.read("epipedobates-boulengei-parcial.fasta","fasta")
print("\n----Descripcion del gen Epipedobates Boulengei----")
print(epipedobates_boulengei)
print("\nSecuencia completa del gen :")
print(epipedobates_boulengei.seq)

epipedobates_espinosai = SeqIO.read("epipedobates-espinosai-parcial.fasta","fasta")
print("\n----Descripcion del gen Epipedobates Espinosai----")
print(epipedobates_espinosai)
print("\nSecuencia completa del gen :")
print(epipedobates_espinosai.seq)


result_handle= NCBIWWW.qblast("blastn", "nt",epipedobates_anthonyi.seq)
print(result_handle)

blast_records = NCBIXML.parse(result_handle)
dir(blast_records)
print(blast_records.__sizeof__)

print("Alineamientos con el nucleotido del Epipedobates anthonyi con la base de datos de BLAST\n")
for b in blast_records:
    for alignment in b.alignments:
        for hsp in alignment.hsps:
            print('\n****Alineamiento****', file=f)
            print('secuencia:', alignment.title, file=f)
            print('longitud:', alignment.length, file=f)
            print('e value:', hsp.expect, file=f)
            print(hsp.query[0:100] + '...', file=f)
            print(hsp.match[0:100] + '...', file=f)
            print(hsp.sbjct[0:100] + '...', file=f)
