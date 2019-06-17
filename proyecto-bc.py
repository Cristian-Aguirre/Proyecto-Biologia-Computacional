
import sys
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtWidgets import QApplication,QMainWindow, QPushButton, QToolTip, QMessageBox, QMenu, QTextEdit
from PyQt5.QtWidgets import QAction, QWidget, QFileDialog, QLabel,QTableWidget,QTableWidgetItem,QVBoxLayout
from PyQt5.QtCore import QCoreApplication, Qt
from PyQt5 import QtGui
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein, generic_rna
from Bio import SeqIO as seq
from Bio import AlignIO as align
import os

"""records=[]
epipedobates_anthonyi = SeqIO.read("epipedobates-anthonyi-parcial.fasta", "fasta")
records.append()
epipedobates_boulengei = SeqIO.read("epipedobates-boulengei-parcial.fasta","fasta")
epipedobates_espinosai = SeqIO.read("epipedobates-espinosai-parcial.fasta","fasta")
epipedobates_darwinwallacei = SeqIO.read("epipedobates-darwinwallacei-parcial.fasta","fasta")
epipedobates_machalilla = SeqIO.read("epipedobates-machalilla-parcial.fasta","fasta")
epipedobates_tricolor = SeqIO.read("epipedobates-tricolor-parcial.fasta","fasta")

seq.write()
"""
records=[]

for filename in os.listdir("secuencias"):
	handle=open("secuencias" + "/"  +filename)
	record=seq.read(handle,"fasta")
	records.append(record)

seq.write(records,"Conjunto_fasta.fasta", "fasta")

aligment=align.read(open("Conjunto_fasta.aln"), "clustal")

print(aligment)

class Prueba(QMainWindow):
    def __init__(self):
        super().__init__()

        self.title="Proyecto Biologia Computacional"
        self.top=0
        self.left=0
        self.width=1200
        self.height=700
        #self.setWindowIcon(QtGui.QIcon(icon))
        button1 = QPushButton("ALINEAR!", self)
        button1.move(300,400)
        button1.setToolTip("<h4>Comprobar sintaxis</h4>")
        button1.clicked.connect(self.onClick1)

        self.InitWindow()


    def InitWindow(self):

        self.label2 = QLabel(self)
        self.label2.setPixmap(QPixmap('images/epipedobates_anthonyi.jpg'))
        self.label2.setGeometry(50,50,250,200)

        nuevaentrada = QAction(QIcon('new_icon.png'), 'Nuevo', self)
        nuevaentrada.setShortcut('Ctrl+E')


        self.toolbar = self.addToolBar('Toolbar')
        self.toolbar.addAction(nuevaentrada)
        self.toolbar.actionTriggered[QAction].connect(self.onClick2)

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
       print("Hola1")

    def onClick2(self):
       print("Hola2")

App = QApplication(sys.argv)
window = Prueba()


sys.exit(App.exec())
