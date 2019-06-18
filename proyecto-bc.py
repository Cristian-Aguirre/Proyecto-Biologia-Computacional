
import sys
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtWidgets import QApplication,QMainWindow, QPushButton, QToolTip, QMessageBox, QMenu, QTextEdit
from PyQt5.QtWidgets import QAction, QWidget, QFileDialog, QLabel,QTableWidget,QTableWidgetItem,QVBoxLayout, QCheckBox
from PyQt5.QtCore import QCoreApplication, Qt
from PyQt5 import QtGui
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein, generic_rna
from Bio import SeqIO as seq
from Bio import AlignIO as align
#from Bio.Blast import NCBIWWW
#from Bio.Blast import NCBIXML
import os
import subprocess
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from Bio import Phylo
import pylab


class Prueba(QMainWindow):
    def __init__(self):
        super().__init__()

        self.title="Proyecto Biologia Computacional"
        self.top=0
        self.left=0
        self.width=1200
        self.height=700
        button1 = QPushButton("ALINEAR!", self)
        button1.move(100,450)
        button1.setToolTip("<h4>Comprobar sintaxis</h4>")
        button1.clicked.connect(self.onClick1)

        button2 = QPushButton("MOSTRAR!", self)
        button2.move(100,400)
        button2.setToolTip("<h4>Agregar secuencias</h4>")
        button2.clicked.connect(self.onClick2)

        button3 = QPushButton("LIMPIAR!", self)
        button3.move(100,500)
        button3.setToolTip("<h4>Limpiar Texto</h4>")
        button3.clicked.connect(self.onClick3)

        button4 = QPushButton("BLAST!", self)
        button4.move(100,550)
        button4.setToolTip("<h4>Alinea con BLAST</h4>")
        button4.clicked.connect(self.onClick4)

        button5 = QPushButton("ARBOL FILOGENETICO!", self)
        button5.setGeometry(50,600,150,30)
        button5.setToolTip("<h4>Muestra el arbol creado</h4>")
        button5.clicked.connect(self.onClick5)

        button6 = QPushButton("EPIPEDOBATES!", self)
        button6.setGeometry(400,10,100,30)
        button6.setToolTip("<h4>Escoger specie</h4>")
        button6.clicked.connect(self.onClick6)

        button7 = QPushButton("COLOSTETHUS!", self)
        button7.setGeometry(700,10,100,30)
        button7.setToolTip("<h4>Escoger familia</h4>")
        button7.clicked.connect(self.onClick7)

        self.InitWindow()

    def InitWindow(self):
        self.contador=0

        self.label2 = QLabel(self)
        self.label2.setPixmap(QPixmap('images/epipedobates_anthonyi.jpg'))
        self.label2.setGeometry(50,30,200,250)
        self.checkbox1 = QCheckBox("Epipedobates\nAnthonyi", self)
        self.checkbox1.setGeometry(50,240,150,50)

        self.label3 = QLabel(self)
        self.label3.setPixmap(QPixmap('images/epipedobates_boulengeri.jpg'))
        self.label3.setGeometry(250,30,450,250)
        self.checkbox2 = QCheckBox("Epipedobates\nBoulengeri",self)
        self.checkbox2.setGeometry(250,220,150,50)

        self.label4 = QLabel(self)
        self.label4.setPixmap(QPixmap('images/epipedobates_darwinwallacei.jpg'))
        self.label4.setGeometry(450,30,650,250)
        self.checkbox3 = QCheckBox("Epipedobates\nDarwinwallacei", self)
        self.checkbox3.setGeometry(450,210,150,50)

        self.label5 = QLabel(self)
        self.label5.setPixmap(QPixmap('images/epipedobates_espinosai.jpg'))
        self.label5.setGeometry(650,30,800,250)
        self.checkbox4 = QCheckBox("Epipedobates\nEspinosai", self)
        self.checkbox4.setGeometry(650,210,150,50)

        self.label6 = QLabel(self)
        self.label6.setPixmap(QPixmap('images/epipedobates_machalilla.jpg'))
        self.label6.setGeometry(850,30,950,250)
        self.checkbox5 = QCheckBox("Epipedobates\nMachalilla",self)
        self.checkbox5.setGeometry(850,210,150,50)

        self.label7 = QLabel(self)
        self.label7.setPixmap(QPixmap('images/epipedobates_tricolor.jpg'))
        self.label7.setGeometry(1050,30,1150,250)
        self.checkbox6 = QCheckBox("Epipedobates\nTricolor", self)
        self.checkbox6.setGeometry(1050,210,150,50)

        self.cajadetexto = QTextEdit(self)
        self.cajadetexto.setGeometry(200,270,1050,430)

        self.setWindowTitle(self.title)
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.blue)
        self.setPalette(p)
        self.setGeometry(self.top,self.left,self.width,self.height)
        self.show()

    def CloseApp(self):
        reply = QMessageBox.question(self,"Mensaje!","¿Seguro de cerrar?",
                QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.close()

    def contextMenuEvent(self, event):
        contextMenu = QMenu(self)
        newAct = contextMenu.addAction("Nuevo")
        openAct = contextMenu.addAction("Abrir")
        quitAct = contextMenu.addAction("Cerrar")

        action = contextMenu.exec_(self.mapToGlobal(event.pos()))

        if action == quitAct:
            self.CloseApp()

    def onClick1(self):
        self.alineamiento_multiple()

    def onClick2(self):
       self.contador=0
       handle_output = open("Conjunto_fasta.fasta","w+")
       if self.checkbox1.isChecked():
            self.f1 = open("secuencias/epipedobates-anthonyi-parcial.fasta","r")
            texto1 = self.f1.read().strip()
            handle_output.write(texto1 +"\n\n")
            self.contador+=1
       if self.checkbox2.isChecked():
            self.f2 = open("secuencias/epipedobates-boulengei-parcial.fasta","r")
            texto2 = self.f2.read().strip()
            handle_output.write(texto2 +"\n\n")
            self.contador+=1
       if self.checkbox3.isChecked():
            self.f3 = open("secuencias/epipedobates-darwinwallacei-parcial.fasta","r")
            texto3 = self.f3.read().strip()
            handle_output.write(texto3 +"\n\n")
            self.contador+=1
       if self.checkbox4.isChecked():
            self.f4 = open("secuencias/epipedobates-espinosai-parcial.fasta","r")
            texto4 = self.f4.read().strip()
            handle_output.write(texto4 +"\n\n")
            self.contador+=1
       if self.checkbox5.isChecked():
            self.f5 = open("secuencias/epipedobates-machalilla-parcial.fasta","r")
            texto5 = self.f5.read().strip()
            handle_output.write(texto5 +"\n\n")
            self.contador+=1
       if self.checkbox6.isChecked():
            self.f6 = open("secuencias/epipedobates-tricolor-parcial.fasta","r")
            texto6 = self.f6.read().strip()
            handle_output.write(texto6 +"\n\n")
            self.contador+=1
       handle_output.close()
       ftotal = open("Conjunto_fasta.fasta","r")
       texto_total = ftotal.read().strip()
       self.cajadetexto.setText(texto_total)

    def alineamiento_multiple(self):
        subprocess.call(['/home/cristian/BC-and-Topicos1/Proyect-BC/procesa_clustal.sh'])
        ft = open("Conjunto_fasta.aln","r")
        input = ft.read().strip()
        self.cajadetexto.setText(input)

    def onClick3(self):
        self.cajadetexto.clear()

    def onClick4(self):
        """print(self.contador)
        if self.checkbox1.isChecked() and self.contador==1:
            output_filename = "output.txt"
            f = open(output_filename,'w')
            sys.oudout = f
            epipedobates_anthonyi = seq.read("secuencias/epipedobates-anthonyi-parcial.fasta", "fasta")
            result_handle= NCBIWWW.qblast("blastn", "nt",epipedobates_anthonyi.seq)
            blast_records = NCBIXML.parse(result_handle)
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
            f1 = open("output.txt","r")
            texto1 = f1.read().strip()
            self.cajadetexto.setText(texto1)
        if self.checkbox2.isChecked() and self.contador==1:
            output_filename = "output.txt"
            f = open(output_filename,'w')
            sys.oudout = f
            epipedobates_boulengei = seq.read("secuencias/epipedobates-boulengei-parcial.fasta", "fasta")
            result_handle= NCBIWWW.qblast("blastn", "nt",epipedobates_boulengei.seq)
            blast_records = NCBIXML.parse(result_handle)
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
            f1 = open("output.txt","r")
            texto1 = f1.read().strip()
            self.cajadetexto.setText(texto1)
"""
    def onClick5(self):
        with open("Conjunto_fasta.aln", "r") as aln:
        #usar AlignIO tpara leer el archivo de alineamiento en formato 'clustal' format
            alignment = align.read(aln, "clustal")
        #calcular la  matriz de distancias
        calculator = DistanceCalculator('identity')
        # añade la matriz de  distancias al objeto calculator y lo retorna
        dm = calculator.get_distance(alignment)

        #Construir el arbol filogenetico aprtir de las distancias
        constructor = DistanceTreeConstructor(calculator)

        upgma_tree = constructor.upgma(dm)

        Phylo.draw_ascii(upgma_tree)
        Phylo.draw(upgma_tree)

    def onClick6(self):
        self.label2.show()
        self.checkbox1.show()
        self.label3.show()
        self.checkbox2.show()
        self.label4.show()
        self.checkbox3.show()
        self.label5.show()
        self.checkbox4.show()
        self.label6.show()
        self.checkbox5.show()
        self.label7.show()
        self.checkbox6.show()

    def onClick7(self):
        self.label2.hide()
        self.checkbox1.hide()
        self.label3.hide()
        self.checkbox2.hide()
        self.label4.hide()
        self.checkbox3.hide()
        self.label5.hide()
        self.checkbox4.hide()
        self.label6.hide()
        self.checkbox5.hide()
        self.label7.hide()
        self.checkbox6.hide()

App = QApplication(sys.argv)
window = Prueba()


sys.exit(App.exec())
