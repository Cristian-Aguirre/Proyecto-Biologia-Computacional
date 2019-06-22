
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
        button6.setGeometry(250,10,100,30)
        button6.setToolTip("<h4>Escoger specie</h4>")
        button6.clicked.connect(self.onClick6)

        button7 = QPushButton("COLOSTETHUS!", self)
        button7.setGeometry(450,10,100,30)
        button7.setToolTip("<h4>Escoger especie</h4>")
        button7.clicked.connect(self.onClick7)

        button8 = QPushButton("AMEEREGA!", self)
        button8.setGeometry(650,10,100,30)
        button8.setToolTip("<h4>Escoger especie</h4>")
        button8.clicked.connect(self.onClick8)

        button9 = QPushButton("PROTEINAS!", self)
        button9.setGeometry(850,10,100,30)
        button9.setToolTip("<h4>Escoger proteinas</h4>")
        button9.clicked.connect(self.onClick9)

        self.InitWindow()

    def InitWindow(self):
        self.contador=0

        self.label2 = QLabel(self)
        self.label2.setPixmap(QPixmap('images/epipedobates_anthonyi.jpg'))
        self.label2.setGeometry(50,30,200,250)
        self.label2.hide()
        self.checkbox1 = QCheckBox("Epipedobates\nAnthonyi", self)
        self.checkbox1.setGeometry(50,240,150,50)
        self.checkbox1.hide()

        self.label3 = QLabel(self)
        self.label3.setPixmap(QPixmap('images/epipedobates_boulengeri.jpg'))
        self.label3.setGeometry(250,30,450,250)
        self.label3.hide()
        self.checkbox2 = QCheckBox("Epipedobates\nBoulengeri",self)
        self.checkbox2.setGeometry(250,220,150,50)
        self.checkbox2.hide()

        self.label4 = QLabel(self)
        self.label4.setPixmap(QPixmap('images/epipedobates_darwinwallacei.jpg'))
        self.label4.setGeometry(450,30,650,250)
        self.label4.hide()
        self.checkbox3 = QCheckBox("Epipedobates\nDarwinwallacei", self)
        self.checkbox3.setGeometry(450,210,150,50)
        self.checkbox3.hide()

        self.label5 = QLabel(self)
        self.label5.setPixmap(QPixmap('images/epipedobates_espinosai.jpg'))
        self.label5.setGeometry(650,30,800,250)
        self.label5.hide()
        self.checkbox4 = QCheckBox("Epipedobates\nEspinosai", self)
        self.checkbox4.setGeometry(650,210,150,50)
        self.checkbox4.hide()

        self.label6 = QLabel(self)
        self.label6.setPixmap(QPixmap('images/epipedobates_machalilla.jpg'))
        self.label6.setGeometry(850,30,950,250)
        self.label6.hide()
        self.checkbox5 = QCheckBox("Epipedobates\nMachalilla",self)
        self.checkbox5.setGeometry(850,210,150,50)
        self.checkbox5.hide()

        self.label7 = QLabel(self)
        self.label7.setPixmap(QPixmap('images/epipedobates_tricolor.jpg'))
        self.label7.setGeometry(1050,30,1150,250)
        self.label7.hide()
        self.checkbox6 = QCheckBox("Epipedobates\nTricolor", self)
        self.checkbox6.setGeometry(1050,210,150,50)
        self.checkbox6.hide()

        self.label8 = QLabel(self)
        self.label8.setPixmap(QPixmap('images/colostethus_fraterdanieli.jpg'))
        self.label8.setGeometry(50,30,200,250)
        self.label8.hide()
        self.checkbox7 = QCheckBox("Colostethus\nFraterdanie", self)
        self.checkbox7.setGeometry(50,220,150,50)
        self.checkbox7.hide()

        self.label9 = QLabel(self)
        self.label9.setPixmap(QPixmap('images/colostethus_inguinalis.jpg'))
        self.label9.setGeometry(310,30,400,250)
        self.label9.hide()
        self.checkbox8 = QCheckBox("Colostethus\nInguinalis", self)
        self.checkbox8.setGeometry(310,210,150,50)
        self.checkbox8.hide()

        self.label10 = QLabel(self)
        self.label10.setPixmap(QPixmap('images/colostethus_agilis.jpg'))
        self.label10.setGeometry(550,30,200,250)
        self.label10.hide()
        self.checkbox9 = QCheckBox("Colostethus\nAgilis", self)
        self.checkbox9.setGeometry(550,210,150,50)
        self.checkbox9.hide()

        self.label11 = QLabel(self)
        self.label11.setPixmap(QPixmap('images/colostethus_panamansis.jpg'))
        self.label11.setGeometry(800,30,200,250)
        self.label11.hide()
        self.checkbox10 = QCheckBox("Colostethus\nPanamansis",self)
        self.checkbox10.setGeometry(800,210,150,50)
        self.checkbox10.hide()

        self.label12 = QLabel(self)
        self.label12.setPixmap(QPixmap('images/colostethus_pratti.jpg'))
        self.label12.setGeometry(1050,30,200,250)
        self.label12.hide()
        self.checkbox11 = QCheckBox("Colostethus\nPratti",self)
        self.checkbox11.setGeometry(1050,210,150,50)
        self.checkbox11.hide()

        self.label13 = QLabel(self)
        self.label13.setPixmap(QPixmap('images/ameerega_bassleri.jpg'))
        self.label13.setGeometry(50,30,200,250)
        self.label13.hide()
        self.checkbox12 = QCheckBox("Ameerega\nBassleri",self)
        self.checkbox12.setGeometry(50,210,150,50)
        self.checkbox12.hide()

        self.label14 = QLabel(self)
        self.label14.setPixmap(QPixmap('images/ameerega_bilinguis.jpg'))
        self.label14.setGeometry(260,30,200,250)
        self.label14.hide()
        self.checkbox13 = QCheckBox("Ameerega\nBilinguis",self)
        self.checkbox13.setGeometry(260,210,150,50)
        self.checkbox13.hide()

        self.label15 = QLabel(self)
        self.label15.setPixmap(QPixmap('images/ameerega_flavopicta.jpg'))
        self.label15.setGeometry(480,30,200,250)
        self.label15.hide()
        self.checkbox14 = QCheckBox("Ameerega\nFlavopicta",self)
        self.checkbox14.setGeometry(480,210,150,50)
        self.checkbox14.hide()

        self.label16 = QLabel(self)
        self.label16.setPixmap(QPixmap('images/ameerega_hahneli.jpg'))
        self.label16.setGeometry(700,55,150,150)
        self.label16.hide()
        self.checkbox15 = QCheckBox("Ameerega\nHahneli",self)
        self.checkbox15.setGeometry(700,210,150,50)
        self.checkbox15.hide()

        self.label17 = QLabel(self)
        self.label17.setPixmap(QPixmap('images/ameerega_macero.jpg'))
        self.label17.setGeometry(880,30,250,250)
        self.label17.hide()
        self.checkbox16 = QCheckBox("Ameerega\nMacero",self)
        self.checkbox16.setGeometry(880,210,150,50)
        self.checkbox16.hide()

        self.label18 = QLabel(self)
        self.label18.setPixmap(QPixmap('images/ameerega_shihuemoy.jpg'))
        self.label18.setGeometry(1100,30,200,250)
        self.label18.hide()
        self.checkbox17 = QCheckBox("Ameerega\nShihuemoy",self)
        self.checkbox17.setGeometry(1100,210,150,50)
        self.checkbox17.hide()


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
       handle_output = open("Conjunto_fasta.fasta","w+")
       if self.checkbox1.isChecked():
            self.f1 = open("secuencias/epipedobates-anthonyi-parcial.fasta","r")
            texto1 = self.f1.read().strip()
            handle_output.write(texto1 +"\n\n")
       if self.checkbox2.isChecked():
            self.f2 = open("secuencias/epipedobates-boulengei-parcial.fasta","r")
            texto2 = self.f2.read().strip()
            handle_output.write(texto2 +"\n\n")
       if self.checkbox3.isChecked():
            self.f3 = open("secuencias/epipedobates-darwinwallacei-parcial.fasta","r")
            texto3 = self.f3.read().strip()
            handle_output.write(texto3 +"\n\n")
       if self.checkbox4.isChecked():
            self.f4 = open("secuencias/epipedobates-espinosai-parcial.fasta","r")
            texto4 = self.f4.read().strip()
            handle_output.write(texto4 +"\n\n")
       if self.checkbox5.isChecked():
            self.f5 = open("secuencias/epipedobates-machalilla-parcial.fasta","r")
            texto5 = self.f5.read().strip()
            handle_output.write(texto5 +"\n\n")
       if self.checkbox6.isChecked():
            self.f6 = open("secuencias/epipedobates-tricolor-parcial.fasta","r")
            texto6 = self.f6.read().strip()
            handle_output.write(texto6 +"\n\n")
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
        print("BLAST en proceso")
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
        self.label8.hide()
        self.checkbox7.hide()
        self.label9.hide()
        self.checkbox8.hide()
        self.label10.hide()
        self.checkbox9.hide()
        self.label11.hide()
        self.checkbox10.hide()
        self.label12.hide()
        self.checkbox11.hide()
        self.label13.hide()
        self.checkbox12.hide()
        self.label14.hide()
        self.checkbox13.hide()
        self.label15.hide()
        self.checkbox14.hide()
        self.label16.hide()
        self.checkbox15.hide()
        self.label17.hide()
        self.checkbox16.hide()
        self.label18.hide()
        self.checkbox17.hide()

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
        self.label8.show()
        self.checkbox7.show()
        self.label9.show()
        self.checkbox8.show()
        self.label10.show()
        self.checkbox9.show()
        self.label11.show()
        self.checkbox10.show()
        self.label12.show()
        self.checkbox11.show()
        self.label13.hide()
        self.checkbox12.hide()
        self.label14.hide()
        self.checkbox13.hide()
        self.label15.hide()
        self.checkbox14.hide()
        self.label16.hide()
        self.checkbox15.hide()
        self.label17.hide()
        self.checkbox16.hide()
        self.label18.hide()
        self.checkbox17.hide()

    def onClick8(self):
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
        self.label8.hide()
        self.checkbox7.hide()
        self.label9.hide()
        self.checkbox8.hide()
        self.label10.hide()
        self.checkbox9.hide()
        self.label11.hide()
        self.checkbox10.hide()
        self.label12.hide()
        self.checkbox11.hide()
        self.label13.show()
        self.checkbox12.show()
        self.label14.show()
        self.checkbox13.show()
        self.label15.show()
        self.checkbox14.show()
        self.label16.show()
        self.checkbox15.show()
        self.label17.show()
        self.checkbox16.show()
        self.label18.show()
        self.checkbox17.show()

    def onClick9(self):
        print("Proteinas")

App = QApplication(sys.argv)
window = Prueba()


sys.exit(App.exec())
