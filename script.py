# coding: latin-1
from astropy.io import ascii

from docutils.nodes import table
from numpy import random
import matplotlib
matplotlib.use("Qt4Agg")

__author__ = 'giovanirech'
import sys
import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import visa
import os
from PyQt4 import QtCore, QtGui
from interface import Ui_MainWindow, _fromUtf8


class StartQT4(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(StartQT4, self).__init__(parent)
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.showMaximized()
        # Conexao de elementos da interface grafica
        self.fig = plt.figure()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        #QtCore.QObject.connect(self.ui.btn_runpoint, QtCore.SIGNAL('clicked()'), self.dynamic_analysis)
        QtCore.QObject.connect(self.ui.btn_detectar, QtCore.SIGNAL('clicked()'), self.instrument_detection)
        QtCore.QObject.connect(self.ui.btn_conectar, QtCore.SIGNAL('clicked()'), self.instrument_connection)
        QtCore.QObject.connect(self.ui.btn_opencal, QtCore.SIGNAL('clicked()'), self.open_calibration)
        QtCore.QObject.connect(self.ui.btn_shortcal, QtCore.SIGNAL('clicked()'), self.short_calibration)
        QtCore.QObject.connect(self.ui.btn_run, QtCore.SIGNAL('clicked()'), self.run_analysis)
        QtCore.QObject.connect(self.ui.toolbtn_savepath, QtCore.SIGNAL('clicked()'), self.get_save_path)
        QtCore.QObject.connect(self.ui.btn_plot, QtCore.SIGNAL('clicked()'), self.plot_test)
        QtCore.QObject.connect(self.ui.btn_plot_permi, QtCore.SIGNAL('clicked()'), self.plot_erei)
        QtCore.QObject.connect(self.ui.btn_plot_RC, QtCore.SIGNAL('clicked()'), self.plot_RC)
        QtCore.QObject.connect(self.ui.btn_plot_ZT, QtCore.SIGNAL('clicked()'), self.plot_ZT)
        QtCore.QObject.connect(self.ui.btn_plot_ZrZi, QtCore.SIGNAL('clicked()'), self.plot_ZrZi)
        QtCore.QObject.connect(self.ui.actionImportar_Dados, QtCore.SIGNAL('triggered()'), self.importar_dados)
        QtCore.QObject.connect(self.ui.blt_skipcal, QtCore.SIGNAL('clicked()'), self.skip_calibration)
        self.is_open_calibrated = False
        self.is_short_calibrated = False
        layout = self.ui.verticalLayout
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)


    # Metodos
    def skip_calibration(self):
        warning = QtGui.QMessageBox()
        warning.setText('Are you sure you want to skip calibration? Without geometry compensation, the collected data'
                        'may be contaminated with the effect of resistivity and capacitance of the accessories, cables and conections.')
        warning.setIcon(QtGui.QMessageBox.Warning)
        warning.addButton('Continue', QtGui.QMessageBox.AcceptRole)
        warning.addButton('Cancel', QtGui.QMessageBox.RejectRole)
        warning.exec_()
        btnclicked = warning.clickedButton().text()
        if btnclicked == 'Continue':
            self.ui.groupBox_analise.setEnabled(True)
            self.ui.groupBox__amostra.setEnabled(True)
            self.ui.btn_run.setEnabled(True)
        else:
            pass

    def importar_dados(self):
        caminho = QtGui.QFileDialog.getOpenFileNameAndFilter(self, 'Select the file', filter='*.txt')
        s = str(caminho[0]).rsplit('/', 1)
        path = s[0]+'/'
        data = np.genfromtxt(str(caminho[0]), skip_header=0)[:, :]
        amostra = 'Identifiction not found'
        diametro = 1
        espessura = 1
        with open(str(caminho[0])) as arquivo:
            for line in arquivo:
                if line.startswith('#D'):
                    diametro = float(line.rstrip('\n').split()[-1])
                if line.startswith('#d'):
                    espessura = float(line.rstrip('\n').split()[-1])
                if line.startswith('#sample'):
                    amostra = ' '.join(line.rstrip('\n').split()[1:-1])
                if line.startswith('#sweep'):
                    if line.rstrip('\n').split()[-1].lower == 'lin':
                        self.ui.var_lin.setChecked(True)
                    else:
                        self.ui.var_log.setChecked(True)
        global Zr, Zi, R, C, er, ei, d, D, Theta, Freq, Z
        nop = len(data[:, 0])
        Freq = data[:, 0]
        Z = data[:, 1]
        Theta = data[:, 2]
        e0 = 8.85418782e-12
        d = espessura*10**(-3)
        D = diametro*10**(-3)
        # er = (2*d*np.cos(np.asarray(Theta)*np.pi/180))/(e0*np.asarray(Z)*np.asarray(Freq)*(D**2)*(np.pi*np.pi)*np.tan(np.asarray(Theta)*np.pi/180))
        A = (np.pi*(D**2))/4
        Zr = Z*np.cos(np.asarray(Theta)*np.pi/180)
        Zi = Z*np.sin(np.asarray(Theta)*np.pi/180)*(-1)
        R = Z/np.cos(np.asarray(Theta)*np.pi/180)
        C = Zi/(Zr*np.asarray(Freq)*R*2*np.pi)

        er = (C*d)/(e0*A)
        ei = er*np.tan((90+np.asarray(Theta))*np.pi/180)
        # ei = er*np.tan((90+np.asarray(Theta))*np.pi/180
        self.ui.btn_plot_ZrZi.setEnabled(True)
        self.ui.btn_plot_ZT.setEnabled(True)
        self.ui.btn_plot_permi.setEnabled(True)
        self.ui.btn_plot_RC.setEnabled(True)
        self.ui.groupBox__amostra.setEnabled(True)
        self.ui.groupBox_analise.setEnabled(True)
        self.ui.amostra_id.setText(amostra)
        self.ui.LineEdit_SavePath.setText(path)
        self.ui.tbox_diametro.setText(str(diametro))
        self.ui.tbox_espessura.setText(str(espessura))
        self.ui.spin_freq_inicial.setValue(int(Freq[0]))
        self.ui.spin_freq_final.setValue(int(Freq[-1]))
        self.ui.num_pontos.setValue(int(nop))
        self.plot_erei()
        self.update_table(np.ndarray.tolist(Freq), np.ndarray.tolist(Z), np.ndarray.tolist(Theta),
                          np.ndarray.tolist(er), np.ndarray.tolist(ei))

    def instrument_detection(self):
        global rm
        is_simulation = self.ui.checkbox_pyvisasim.checkState()
        if is_simulation:
            rm = visa.ResourceManager('@sim')
        else:
            rm = visa.ResourceManager()
        equipaments = list(rm.list_resources())
        self.ui.combobox_equipamentos.addItems(equipaments)

    def instrument_connection(self):
        global instrument
        address = self.ui.combobox_equipamentos.currentText()
        try:
            instrument = rm.open_resource(unicode(address))
        except:
            message = QtGui.QMessageBox()
            message.setText('Conection error.')
            message.setDetailedText('Please, verify that the equipament is turned on and connected to the computer.'
                                     'Make sure that the conection cable is working and you have installed all the necessary drivers.')
            message.exec_()
            return
        instrument_name = instrument.query('*IDN?')
        self.ui.label_inst_name.setText(instrument_name)
        self.ui.group_calibration.setEnabled(True)


    def open_calibration(self):
        instrument.write('HOLD')
        warning = QtGui.QMessageBox()
        warning.setText('Please, set the accessory to the OPEN configuration')
        warning.setIcon(QtGui.QMessageBox.Information)
        warning.addButton('Continue', QtGui.QMessageBox.AcceptRole)
        warning.addButton('Cancel', QtGui.QMessageBox.RejectRole)
        warning.exec_()
        response = warning.clickedButton().text()
        if response == 'Continue':
            instrument.write('E4TP OFF')
            instrument.write('TRGS INT')
            instrument.write('ESNB 1')
            instrument.write('*SRE 4')
            instrument.write('*CLS')
            instrument.write('COMA')
            del instrument.timeout
            instrument.write('*WAI')
            instrument.ask('*OPC?')
            complete = QtGui.QMessageBox()
            complete.setText('OPEN compensation complete.')
            complete.setIcon(QtGui.QMessageBox.Information)
            complete.addButton('OK', QtGui.QMessageBox.AcceptRole)
            complete.exec_()
            self.is_open_calibrated = True
            if self.is_open_calibrated and self.is_short_calibrated:
                self.ui.groupBox_analise.setEnabled(True)
                self.ui.groupBox__amostra.setEnabled(True)
                self.ui.btn_run.setEnabled(True)
        else:
            pass

    def short_calibration(self):
        instrument.write('HOLD')
        warning = QtGui.QMessageBox()
        warning.setText('Please, set the accessory to the SHORT configuration')
        warning.setIcon(QtGui.QMessageBox.Information)
        warning.addButton('Continue', QtGui.QMessageBox.AcceptRole)
        warning.addButton('Cancel', QtGui.QMessageBox.RejectRole)
        warning.exec_()
        response = warning.clickedButton().text()
        if response == 'Continue':
            instrument.write('E4TP OFF')
            instrument.write('TRGS INT')
            instrument.write('ESNB 1')
            instrument.write('*SRE 4')
            instrument.write('*CLS')
            instrument.write('COMB')
            del instrument.timeout
            instrument.write('*WAI')
            instrument.ask('*OPC?')
            complete = QtGui.QMessageBox()
            complete.setText('SHORT compensation complete.')
            complete.setIcon(QtGui.QMessageBox.Information)
            complete.addButton('OK', QtGui.QMessageBox.AcceptRole)
            complete.exec_()
            self.is_short_calibrated = True
            if self.is_open_calibrated and self.is_short_calibrated:
                self.ui.groupBox_analise.setEnabled(True)
                self.ui.groupBox__amostra.setEnabled(True)
                self.ui.btn_run.setEnabled(True)
        else:
            pass

    def get_save_path(self):
        path = QtGui.QFileDialog.getExistingDirectory(self, 'Select destination folder')
        self.ui.LineEdit_SavePath.setText(path)

    def no_blank_fields(self):
        if str(self.ui.amostra_id.text()) == '':
            warning = QtGui.QMessageBox()
            warning.setIcon(QtGui.QMessageBox.Warning)
            warning.setText('You need t provide a identification to your sample.')
            warning.addButton('Ok', QtGui.QMessageBox.AcceptRole)
            warning.exec_()
            self.ui.amostra_id.hasFocus()
            return False
        if str(self.ui.LineEdit_SavePath.text()) == '':
            warning = QtGui.QMessageBox()
            warning.setIcon(QtGui.QMessageBox.Warning)
            warning.setText('Select a destination folder to save the colected data.')
            warning.addButton('Ok', QtGui.QMessageBox.AcceptRole)
            warning.exec_()
            return False
        if str(self.ui.tbox_diametro.text()) == '' or str(self.ui.tbox_espessura.text()) == '':
            warning = QtGui.QMessageBox()
            warning.setIcon(QtGui.QMessageBox.Warning)
            warning.setText('You need to provide the sample dimensions for permittivity and capacitance computations')
            warning.addButton('Ok', QtGui.QMessageBox.AcceptRole)
            warning.exec_()
            self.ui.tbox_diametro.hasFocus()
            return False
        return True


    def run_analysis(self):
        try:
            instrument
        except NameError:
            errorbox = QtGui.QMessageBox()
            errorbox.setIcon(QtGui.QMessageBox.Critical)
            errorbox.setText('Instrument no recognized or not connected.')
            errorbox.addButton('Ok', QtGui.QMessageBox.AcceptRole)
            errorbox.exec_()
        if not(self.no_blank_fields()):
            return
        if self.ui.var_lin.isChecked():
            sweepcmd = 'LIN'
        else:
            sweepcmd = 'LOG'

        if self.ui.checkbox_ptavg.isChecked():
            paver = 'ON'
        else:
            paver = 'OFF'
        id_amostra = self.ui.amostra_id.text()
        diametro = float(self.ui.tbox_diametro.text())
        espessura = float(self.ui.tbox_espessura.text())
        instrument.write('STAR %s' %self.ui.spin_freq_inicial.value())
        instrument.write('STOP %s' %self.ui.spin_freq_final.value())
        instrument.write('POIN %s' %self.ui.num_pontos.value())
        instrument.write('POWMOD VOLT')
        instrument.write('POWE %s' %self.ui.spinbox_tensao.value())
        instrument.write('SWPT %s' %sweepcmd)
        instrument.write('BWFACT %s' %self.ui.spinbox_banda.value())
        instrument.write('PAVERFACT %s' %self.ui.spin_ptavg.value())
        instrument.write('PAVER %s' %paver)
        instrument.write('MEASTAT %s' %paver)
        instrument.write('MEAS IMPH')
        instrument.write('SING')
        instrument.write('TRGS INT')
        instrument.write('ESNB 1')
        instrument.write('*SRE 4')
        instrument.write('*CLS')
        del instrument.timeout
        instrument.write('*WAI')
        instrument.ask('*OPC?')
        instrument.write('TRAC A')
        instrument.write('AUTO')
        instrument.write('TRAC B')
        instrument.write('AUTO')
        global Zr, Zi, R, C, er, ei, d, D, Theta, Freq, Z
        instrument.write('FORM5')
        instrument.write('TRAC A')
        nop = instrument.ask('POIN?')
        Z_data = instrument.ask_for_values('OUTPDTRC?')
        Z = range(int(nop))
        for x in range(0, int(nop)):
            Z[x] = Z_data[2*x]
        instrument.write('TRAC B')
        Theta_data = instrument.ask_for_values('OUTPDTRC?')
        Theta = range(int(nop))
        for x in range(0, int(nop)):
            Theta[x] = Theta_data[2*x]
        Freq = instrument.ask_for_values('OUTPSWPRM?')
        #stats = instrument.ask_for_values('MEASTAT?')
        e0 = 8.85418782e-12
        d = espessura*10**(-3)
        D = diametro*10**(-3)
        # er = (2*d*np.cos(np.asarray(Theta)*np.pi/180))/(e0*np.asarray(Z)*np.asarray(Freq)*(D**2)*(np.pi*np.pi)*np.tan(np.asarray(Theta)*np.pi/180))
        A = (np.pi*(D**2))/4
        Zr = Z*np.cos(np.asarray(Theta)*np.pi/180)
        Zi = Z*np.sin(np.asarray(Theta)*np.pi/180)*(-1)
        R = Z/np.cos(np.asarray(Theta)*np.pi/180)
        C = Zi/(Zr*np.asarray(Freq)*R*2*np.pi)

        er = (C*d)/(e0*A)
        ei = er*np.tan((90+np.asarray(Theta))*np.pi/180)
        # ei = er*np.tan((90+np.asarray(Theta))*np.pi/180
        self.plot_erei()
        self.update_table(Freq, Z, Theta, np.ndarray.tolist(er), np.ndarray.tolist(ei))
        path = str(self.ui.LineEdit_SavePath.text())
        filename = str(id_amostra)
        M = np.c_[np.asarray(Freq), np.asarray(Z), np.asarray(Theta), er, ei]
        np.savetxt('%s\%s.txt' % (path, filename), M, fmt='%1.4e',
                   delimiter='\t',
                   header='#Freq(Hz)\t Z(ohms)\t Phase(degrees) \t er_Re \t er_Im',
                   comments='#Sample %s\n#d(mm)  %s\n#D(mm) %s \n#Sweep  %s \n#Voltage  %s'
                            % (id_amostra, espessura, diametro, sweepcmd, self.ui.spinbox_tensao.value()))
        self.ui.btn_plot_ZrZi.setEnabled(True)
        self.ui.btn_plot_ZT.setEnabled(True)
        self.ui.btn_plot_permi.setEnabled(True)
        self.ui.btn_plot_RC.setEnabled(True)

    def dynamic_analysis(self):
        #message = QtGui.QMessageBox()
        #message.setText('Fun\E7\E3o ainda em constru\E7\E3o.')
        #message.setDetailedText('Por favor, utilize o modo padr\E3o de an\E1lise por enquanto.')
        #message.exec_()
        #return
        try:
            instrument
        except NameError:
            errorbox = QtGui.QMessageBox()
            errorbox.setIcon(QtGui.QMessageBox.Critical)
            errorbox.setText('Instrument not connected or not recognized.')
            errorbox.addButton('Ok', QtGui.QMessageBox.AcceptRole)
            errorbox.exec_()
        if not(self.no_blank_fields()):
            return
        if self.ui.var_lin.isChecked():
            sweepcmd = 'LIN'
        else:
            sweepcmd = 'LOG'

        if self.ui.checkbox_ptavg.isChecked():
            paver = 'ON'
        else:
            paver = 'OFF'
        id_amostra = self.ui.amostra_id.text()
        diametro = float(self.ui.tbox_diametro.text())
        espessura = float(self.ui.tbox_espessura.text())
        instrument.write('STAR %s' %self.ui.spin_freq_inicial.value())
        instrument.write('STOP %s' %self.ui.spin_freq_final.value())
        instrument.write('POIN %s' %self.ui.num_pontos.value())
        instrument.write('POWMOD VOLT')
        instrument.write('POWE %s' %self.ui.spinbox_tensao.value())
        instrument.write('SWPT %s' %sweepcmd)
        instrument.write('BWFACT %s' %self.ui.spinbox_banda.value())
        instrument.write('PAVERFACT %s' %self.ui.spin_ptavg.value())
        instrument.write('PAVER %s' %paver)
        instrument.write('MEASTAT %s' %paver)
        instrument.write('MEAS IMPH')
        instrument.write('MAN')
        instrument.write('TRGS INT')
        instrument.write('ESNB 1')
        instrument.write('*SRE 4')
        instrument.write('*CLS')
        del instrument.timeout
        instrument.write('*WAI')
        instrument.ask('*OPC?')
        instrument.write('TRAC A')
        instrument.write('AUTO')
        instrument.write('TRAC B')
        instrument.write('AUTO')
        global Zr, Zi, R, C, er, ei, d, D, Theta, Freq, Z
        instrument.write('FORM5')
        instrument.write('TRAC A')
        nop = instrument.ask('POIN?')
        Z_data = instrument.ask_for_values('OUTPDTRC?')
        Z = range(int(nop))
        for x in range(0, int(nop)):
            Z[x] = Z_data[2*x]
        instrument.write('TRAC B')
        Theta_data = instrument.ask_for_values('OUTPDTRC?')
        Theta = range(int(nop))
        for x in range(0, int(nop)):
            Theta[x] = Theta_data[2*x]
        Freq = instrument.ask_for_values('OUTPSWPRM?')
        #stats = instrument.ask_for_values('MEASTAT?')
        e0 = 8.85418782e-12
        d = espessura*10**(-3)
        D = diametro*10**(-3)
        # er = (2*d*np.cos(np.asarray(Theta)*np.pi/180))/(e0*np.asarray(Z)*np.asarray(Freq)*(D**2)*(np.pi*np.pi)*np.tan(np.asarray(Theta)*np.pi/180))
        A = (np.pi*(D**2))/4
        Zr = Z*np.cos(np.asarray(Theta)*np.pi/180)
        Zi = Z*np.sin(np.asarray(Theta)*np.pi/180)*(-1)
        R = Z/np.cos(np.asarray(Theta)*np.pi/180)
        C = Zi/(Zr*np.asarray(Freq)*R*2*np.pi)

        er = (C*d)/(e0*A)
        ei = er*np.tan((90+np.asarray(Theta))*np.pi/180)
        # ei = er*np.tan((90+np.asarray(Theta))*np.pi/180
        self.plot_erei()
        self.update_table(Freq, Z, Theta, np.ndarray.tolist(er), np.ndarray.tolist(ei))
        path = str(self.ui.LineEdit_SavePath.text())
        filename = str(id_amostra)
        M = np.c_[np.asarray(Freq), np.asarray(Z), np.asarray(Theta), er, ei]
        np.savetxt('%s\%s.txt' % (path, filename), M, fmt='%1.4e',
                   delimiter='\t',
                   header='#Freq(Hz)\t Z(ohms)\t Phase(degrees) \t er_Re \t er_Im',
                   comments='#sample %s\n#d(mm) = %s\n#D(mm) = %s \n#sweep = %s \n#voltage = %s'
                            % (id_amostra, espessura, diametro, sweepcmd, self.ui.spinbox_tensao.value()))
        self.ui.btn_plot_ZrZi.setEnabled(True)
        self.ui.btn_plot_ZT.setEnabled(True)
        self.ui.btn_plot_permi.setEnabled(True)
        self.ui.btn_plot_RC.setEnabled(True)

        # path = str(self.ui.LineEdit_SavePath.text())
        # filename = str(id_amostra)
        # M = np.c_[np.asarray(Freq), np.asarray(Z), np.asarray(Theta), er, ei]
        # np.savetxt('%s\%s.txt' % (path, filename), M, fmt='%1.4e',
        #           delimiter='\t',
        #           header='Freq(Hz)\t Z(ohms)\t Phase(degrees) \t er_Re \t er_Im',
        #           comments='#Amostra %s\n #d = %s \n #D = %s \n' % (id_amostra, espessura, diametro))

    def plot_test(self):
        data = [random.random() for i in range(10)]

        # create an axis
        ax = self.fig.add_subplot(111)

        # discards the old graph
        ax.hold(False)

        # plot data
        ax.plot(data, '*-')

        # refresh canvas
        self.canvas.draw()
        self.update_table(data, data, data, data, data)

    def plot_erei(self):
        self.plot_data(Freq, 'Frequency (Hz)', er, '$\epsilon\'$', ei, '$\epsilon\'\'$', 'log')

    def plot_RC(self):
        self.plot_data(Freq, 'Frequency (Hz)', R, 'R (Ohms)', C, 'C (F)', 'log')

    def plot_ZT(self):
        self.plot_data(Freq, 'Frequency (Hz)', Z, 'Impedance (Ohms)', Theta, 'Phase (degrees)', 'log')

    def plot_ZrZi(self):
        self.fig.clear()
        self.ax1 = self.fig.add_subplot(111)
        self.ax1.hold(False)
        self.ax1.plot(Zr, Zi, 'blue')
        self.ax1.set_title(self.ui.amostra_id.text(), fontsize='12')
        self.ax1.grid(which='both')
        self.ax1.set_xlabel('Real Impedance', fontsize='12')
        self.ax1.set_ylabel('Imaginary Impedance', fontsize='18', color='blue')
        plt.xlim(min(Zr), max(Zr))
        self.canvas.draw()

    def plot_data(self, x, labelx, y1, labely1, y2, labely2, xscale):
        """

        :rtype : object
        """
        self.fig.hold(False)
        self.fig.clear()
        self.ax1 = self.fig.add_subplot(111)
        self.ax1.hold(False)
        self.ax1.plot(x, y1, 'blue')
        self.ax1.set_title(self.ui.amostra_id.text(), fontsize='12')
        self.ax2 = self.ax1.twinx()
        self.ax2.plot(x, y2, 'green')
        self.ax1.grid(which='both')
        self.ax1.set_xlabel(labelx, fontsize='12')
        self.ax1.set_ylabel(labely1, fontsize='18', color='blue')
        self.ax2.set_ylabel(labely2, fontsize='18', color='green')
        plt.xlim(x[0], x[-1])
        self.ax1.set_xscale(xscale)
        self.ax2.set_xscale(xscale)
        for i in self.ax1.get_yticklabels():
            i.set_color('blue')
        for i in self.ax2.get_yticklabels():
            i.set_color('green')
        # plt.savefig('%s.png' %file_name.value, dpi = 300)
        # self.ui.plot_window.
        self.canvas.draw()

    def update_table(self, freq, imped, teta, er, ei):
        data = [freq, imped, teta, er, ei]
        header = ['Frequency (Hz)', 'Impedance (Ohms)', 'Phase (deg)', 'Real Permittivity', 'Imaginary Permittivity']
        self.table_model = MyTableModel(data, header, parent=self)
        self.ui.tableView.setModel(self.table_model)
        # set row height
        nrows = len(data[0])
        for row in xrange(nrows):
            self.ui.tableView.setRowHeight(row, 18)
        for col in [1, 3, 4]:
            self.ui.tableView.resizeColumnToContents(col)

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, datain, headerdata, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = datain
        self.headerdata = headerdata
        self.rows = range(0, len(datain[0]))

    def rowCount(self, parent):
        return len(self.arraydata[0])

    def columnCount(self, parent):
        return len(self.arraydata)

    def data(self, index, role):
        if not index.isValid():
            return QtCore.QVariant()
        elif role != QtCore.Qt.DisplayRole:
            return QtCore.QVariant()
        return QtCore.QVariant(self.arraydata[index.column()][index.row()])

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(self.headerdata[col])
        return QtCore.QVariant()

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = StartQT4()
    myapp.show()
    sys.exit(app.exec_())
