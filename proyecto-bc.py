
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
import os
import subprocess


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

        button4 = QPushButton("ARBOL FILOGENETICO!", self)
        button4.setGeometry(50,550,150,30)
        button4.setToolTip("<h4>Muestra el arbol creado</h4>")
        button4.clicked.connect(self.onClick4)

        self.InitWindow()

    def InitWindow(self):

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
       handle_output = open("Conjunto_fasta.fasta","w+")
       if self.checkbox1.isChecked():
            f1 = open("secuencias/epipedobates-anthonyi-parcial.fasta","r")
            texto1 = f1.read().strip()
            handle_output.write(texto1 +"\n\n")
       if self.checkbox2.isChecked():
            f2 = open("secuencias/epipedobates-boulengei-parcial.fasta","r")
            texto2 = f2.read().strip()
            handle_output.write(texto2 +"\n\n")
       if self.checkbox3.isChecked():
            f3 = open("secuencias/epipedobates-darwinwallacei-parcial.fasta","r")
            texto3 = f3.read().strip()
            handle_output.write(texto3 +"\n\n")
       if self.checkbox4.isChecked():
            f4 = open("secuencias/epipedobates-espinosai-parcial.fasta","r")
            texto4 = f4.read().strip()
            handle_output.write(texto4 +"\n\n")
       if self.checkbox5.isChecked():
            f5 = open("secuencias/epipedobates-machalilla-parcial.fasta","r")
            texto5 = f5.read().strip()
            handle_output.write(texto5 +"\n\n")
       if self.checkbox6.isChecked():
            f6 = open("secuencias/epipedobates-tricolor-parcial.fasta","r")
            texto6 = f6.read().strip()
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
        print("HOLA")

App = QApplication(sys.argv)
window = Prueba()


sys.exit(App.exec())
