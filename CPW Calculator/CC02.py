'''
Program to calculate the basic characteristics of Coplanar Waveguide Resonators, CPW.
Using Geometrical factors describing the CPW, the superconductor used and the
dielectric loss the program calculate the Quality factor and the resonance
frequency.

Author: Francisco Rouxinol - 07/2014
'''


# Python folder 
#!/opt/local/bin/python2.7


import sys
from math import *
from numpy import *
from scipy.special import ellipk
from PyQt4 import QtCore, QtGui
from window import Ui_MainWindow
from ui_main import Ui_Form

#Constants

cl=299792458 # speed light
e0=8.8541878e-12 # 
u0=1.25664E-06 # 
e1=11.6 # Silicon dielectric constant at 1K
E0=1 # vaccum 
#Tg=5.00E-08 # (T=100mK)	
Rl=50 # load resistence
h=6.602E-34 #Planck
Ccap=logspace(-18,-13,30)
#################################

#calculation

def cal(w,s,t,h,l,e,rho,tc,ck,Tg):
    
    #__________Dielectric silicon air interface_____________________ 
    
    k0  = w/(w+2*s)
    k0l =sqrt(1-k0**2)
    k1  = (sinh(pi*w/(4*h)))/(sinh(pi*(2*s+w)/(4*h)))
    k1l = sqrt(1-k1**2)
    
    K0  = ellipk(k0)
    K0l = ellipk(k0l)
    K1  = ellipk(k1)
    K1l = ellipk(k1l)
    
    eeff=1+(e-1)*(K1*K0l)/(2*K1l*K0)
    
    #__________Kinetic inductance_______________________________
    
    a = (1/(2*k0*k0*K0*K0))
    b = -log(t/(4*w))
    c = -k0*log(t/(4*(w+2*s)))
    d = (2*(w+s)/(w+2*s))*log(s/(w+s))
    g = a*(b+c+d)
    #print g, a, b, c , d
    l0 = 1.05e-3*sqrt(rho/tc)
    Lk = u0*g*(l0**2)/(w*t)
    #___________Circuit values________________________________________
    
    Ll = (u0*K0l/(K0*4))+Lk
    Cl =4*e0*eeff*K0/K0l
    f0=cl/(sqrt(eeff)*2*l)
    z0=sqrt(Ll/Cl)
    #print " a", e0, 'b', eeff, 'c', K0,'d', K0l
    
    #____________Dielectric Loss_____________________________________
    
    lambda0 = 2*pi*f0*sqrt(eeff)/cl
    alfa=(e1/sqrt(eeff))*((eeff-1)/(e1-1))*pi*Tg*lambda0/2
    #print alfa
    # __________Circuit Values with loss___________________________
    #print l
    L = 2*Ll*l/(pi*pi)
    C = Cl*l/2
    #print Cl, C
    R = z0/(alfa*l)
    F0 = 1/(2*pi*sqrt(L*C))
    wn = 2*pi*F0
    Qi = 2*pi*F0*R*C

    #print wn, 2*pi*F0
    #______________Final Device parameters____________________________
    
    a = (wn**2)*(ck**2)*Rl
    Rstar = (1+a*Rl)/a
    Cstar = (ck)/(1+a*Rl)
    Wstar = 1/(sqrt(L*(C+2*ck)))
    fstar = Wstar/(2*pi)

   # Qradiation = 5e-3* (cl/(fstar*w))**2
    #Q1a = Qradiation*Qi/(Qi+Qradiation)

    Qext = wn*Rstar*C/2
    #Qload = Qext*Qi/(Qext+Qi)
    Qload = wn*(C+2*Cstar)/(1/R+2/Rstar)
    
    
    return eeff,l0,Lk,Ll,Cl,f0,alfa, z0,L,C,R,F0,Qi,Rstar,Cstar,fstar,Qext,Qload
# Test plot
#from ui_main import Ui_Form
#
#class Plot_Widget(QWidget,  Ui_Form):
#    def __init__(self, data2plot=None, parent = None):
#        super(Plot_Widget, self).__init__(parent)
#        self.setupUi(self)
#        
#        QObject.connect(self.plotBtn, SIGNAL("clicked()"),self.plotData)
#        
#        
#    def plotData(self):
#        x, y = S.rand(2, 30)
#        self.plotWidget.canvas.ax.plot(x, y, 'o')
#        self.plotWidget.canvas.draw()

####TEST
#
#class Plot_Widget(QtGui.QMainWindow, Ui_Form):
#    def __init__(self, data2plot=None, parent = None):
#        super(Plot_Widget, self).__init__(parent)
#        self.setupUi(self)
#        
#        QtCore.QObject.connect(self.plotBtn, QtCore.SIGNAL("clicked()"),self.plotData)
#        
#         
#    def plotData(self):
#        x, y = rand(2, 30)
#        self.plotWidget.canvas.ax.plot(x, y, 'o')
#        self.plotWidget.canvas.draw()
#### END test



# open window
class StartQT4(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        
        #connect to run button
        
        QtCore.QObject.connect(self.ui.runButton, QtCore.SIGNAL(("clicked()")), self.handleCalculate)
        QtCore.QObject.connect(self.ui.plotBtn,   QtCore.SIGNAL(("clicked()")), self.plotData)
      #  QtCore.QObject.connect(self.ui.plotBtn,   QtCore.SIGNAL(("clicked()")), self.handleCalculate)
       
       
#    def plotData(self):
#        x, y = rand(2, 30)
#        self.ui.plotWidget.canvas.ax.plot(x, y, 'o')
#        self.ui.plotWidget.canvas.draw()
 
    def plotData(self):
        self.handleCalculate()
        w   = float(str(self.ui.wEdit.text()))
        s   = float(str(self.ui.sEdit.text()))
        t   = float(str(self.ui.tEdit.text()))
        h   = float(str(self.ui.hEdit.text()))
        l   = float(str(self.ui.lEdit.text()))
        e   = float(str(self.ui.eEdit.text()))
        rho = float(str(self.ui.rhoEdit.text()))
        tc  =  float(str(self.ui.tcEdit.text()))
        ck  = float(str(self.ui.ccEdit.text()))
        Tg  = float(str(self.ui.tgEdit.text()))
        ck = []
        Ql = []
        for n in range(0,len(Ccap)):
           # print 'n = ', n
           # print 'Ccap = ', Ccap[n]
            
            ck.append(Ccap[n])
           # print 'ck=  ', ck[n] eeff,l0,Lk,Ll,Cl,f0,alfa, z0,L,C,R,F0,Qi,Rstar,Cstar,fstar,Qext,
            
            Qload = cal(w,s,t,h,l,e,rho,tc,ck[n],Tg)[17] 
            Ql.append(Qload)
            
        savetxt('Ql.txt',Ql)
            #print 'end'
        
        self.ui.plotWidget.canvas.ax.plot(ck, Ql, marker='o', linestyle='-', color='r')
#            axes[0].plot(ng_vec, energies[:,n])
            #self.ui.plotWidget.canvas.ax.set_ylim(-10, ymax[0])
        self.ui.plotWidget.canvas.ax.set_yscale('log')
        self.ui.plotWidget.canvas.ax.set_xscale('log')
        self.ui.plotWidget.canvas.ax.set_xlabel(r'$C_c$', fontsize=12)
        self.ui.plotWidget.canvas.ax.set_ylabel(r'$Q_l$', fontsize=12)
            
            #self.ui.plotWidget.canvas.ax.plot(x, y, 'o')
       #     self.ui.plotWidget.canvas.draw()
        self.ui.plotWidget.canvas.draw()

### TEst
# #       def newWindow(self):
#        self.myOtherWindow = Plot_Widget()
#        self.myOtherWindow.show()
#        #app = QApplication(sys.argv)
#        #plot = Plot_Widget()#data2plot,  varDict = totaldict)
#        #plot.show()
#         #sys.exit(app.exec_()) 
### TEST


            
    def handleCalculate(self):
        w   = float(str(self.ui.wEdit.text()))
        s   = float(str(self.ui.sEdit.text()))
        t   = float(str(self.ui.tEdit.text()))
        h   = float(str(self.ui.hEdit.text()))
        l   = float(str(self.ui.lEdit.text()))
        e   = float(str(self.ui.eEdit.text()))
        rho = float(str(self.ui.rhoEdit.text()))
        tc  =  float(str(self.ui.tcEdit.text()))
        ck  = float(str(self.ui.ccEdit.text()))
        Tg  = float(str(self.ui.tgEdit.text()))
        #print ck
        
        
        
        
        eeff,l0,Lk,Ll,Cl,f0,alfa, z0,L,C,R,F0,Qi,Rstar,Cstar,fstar,Qext,Qload = cal(w,s,t,h,l,e,rho,tc,ck,Tg)
        #print ans
        self.ui.eeffLabel.setText( str("%.3f" % (eeff)))
        self.ui.lambdaLabel.setText( str("%.3E" % (l0)))
        self.ui.kineticLabel.setText('%.3E H/m' % (Lk)) 
        self.ui.alfaLabel.setText( str("%.3E" % (alfa)))
        
        self.ui.fEdit.setText(str("%.3E" % (f0)))
        self.ui.llEdit.setText(str("%.3E" % ((Ll))))
        self.ui.clEdit.setText(str("%.3E" % ((Cl))))
        self.ui.z0Edit.setText(str("%.2f" % ((z0))))
        
        self.ui.LEdit.setText(str("%.2E" % ((L))))
        self.ui.CEdit.setText(str("%.2E" % ((C))))
        self.ui.REdit.setText(str("%.2E" % ((R))))
        self.ui.f0Edit.setText(str("%.2E" % ((F0))))
        self.ui.qiEdit.setText(str("%.2f" % ((Qi/1e3))))
        self.ui.qeEdit.setText(str("%.2f" % ((Qext/1e3))))
        self.ui.qlEdit.setText(str("%.2f" % ((Qload/1e3))))
        self.ui.frEdit.setText(str("%.3E" % ((fstar))))
        
        #self.ui.kineticLabel.setText('L Kinetic = ' +str("%.3E" % (Lk)))
        #self.ui.fEdit.setText(str(ans))
        

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = StartQT4()
    myapp.show()
    sys.exit(app.exec_())
    
    
    
    