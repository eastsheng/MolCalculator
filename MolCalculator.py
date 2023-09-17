"""
# 打包
pyinstaller -D -w ./MolCalculator.py --clean -i ./imgs/icons-64.png
# 一个分子计算器，可以实现：
# 1. 输入化学式，输出摩尔质量
# 2. 输入分子化学式和分子数目（至少两组），输出质量分数
# 3. 输入分子化学式和分子数目（至少两组），输出摩尔分数
# 4. 输入分子化学式和分子数目，与体系大小(Å)，输出质量密度g/mL
# 5. 添加Draw molecules, MolCalc和Molview三个WEB工具

"""
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QHBoxLayout ,QVBoxLayout,\
QTextEdit, QMenuBar, QMenu, QAction, QMessageBox, QWidget, QVBoxLayout, QLabel, QDialog
from PyQt5.QtGui import QColor, QFont, QIcon, QPixmap
from PyQt5.QtCore import Qt, QCoreApplication, QSize, QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView
import scripts.cal_fraction as cf


px, py = 100,100
x, y = 600, 350
mx, my = 1920, 1080
m = 10

class AboutDialog(QDialog):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("MolCalculator")
        self.setWindowIcon(QIcon("./imgs/icons-64.png"))
        self.setGeometry(px+100, py+100, int(x/2), int((y+m)/3))
        self.setMinimumSize(int(x/2), int((y+m)/3))
        layout = QVBoxLayout()

        title = QLabel("MolCalculator")
        title.setStyleSheet("font-size: 16px; color: #333333;font-family: Arial;")
        title.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
        layout.addWidget(title)

        description = QLabel()
        description.setOpenExternalLinks(True)
        description.setText('A molecular calculator<br>Copyright @ 2023 MolCalculator<br>written by <a href="https://github.com/eastsheng/MolCalculator">eastsheng</a>')
        description.setStyleSheet("font-size: 12px; color: #333333;font-family: Arial;")
        description.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
        layout.addWidget(description)

        self.setLayout(layout)
        self.setStyleSheet("""
            QDialog {
                background-color: #f5f5f5;  /* 设置背景颜色 */
                border: 1px solid #f2f2f2;  /* 设置边框 */
            }
            QLabel {
                margin: 2px;
            }
        """)


class ChemicalCalculator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):        
        self.setWindowTitle("MolCalculator")
        self.setWindowIcon(QIcon("./imgs/icons-64.png"))
        self.setGeometry(px, py, x, y+m)
        # self.setMinimumSize(x, y+m)
        # self.setMaximumSize(mx, my)
        self.setStyleSheet("background-color: #f0f0f0;")

        self.label = QtWidgets.QLabel(self)
        self.label.setText("Formula:")
        self.label.move(20, 30+m)
        self.label.setStyleSheet("font-size: 20px;font-family: Arial;")

        self.textbox = QtWidgets.QLineEdit(self)
        self.textbox.setText("NaCl")
        self.textbox.setGeometry(130, 30+m, 200, 30)
        self.textbox.setStyleSheet("background-color: white; font-size: 20px;border: 1px solid #ccc;")
        self.textbox.setStyleSheet("""
            QLineEdit {
                background-color: #f5f5f5;  /* 设置背景颜色 */
                font-family: Arial;  /* 设置字体样式 */
                font-size: 16px;
                border: 1px solid #ccc;  /* 设置边框 */
            }
        """)
        self.button = QtWidgets.QPushButton(self)
        self.button.setText("Run")
        self.button.setGeometry(350, 30+m, 55, 30)
        self.button.setStyleSheet("font-family: Arial; background-color: #808080; color: white; font-size: 16px;")
        self.button.clicked.connect(self.calculate_molecular_weight)
        
        self.out_textbox = QTextEdit(self)
        self.out_textbox.setGeometry(20, 70+m, 385, 220)
        self.out_textbox.setReadOnly(True)
        self.out_textbox.setStyleSheet("""
            QTextEdit {
                background-color: #f5f5f5;  /* 设置背景颜色 */
                font-family: Arial;  /* 设置字体样式 */
                font-size: 20px;
                border: 1px solid #ccc;  /* 设置边框 */
            }
        """)
        self.setStyleSheet("QPushButton:hover{background-color: #808080;}")

        self.button_clean = QtWidgets.QPushButton(self)
        self.button_clean.setText("Clean")
        self.button_clean.setGeometry(60, y-40+m, 60, 30)
        self.button_clean.setStyleSheet("font-family: Arial;background-color: #495c69; color: white; font-size: 16px;")
        self.button_clean.clicked.connect(self.clear_output_textbox)

        self.button_exit = QtWidgets.QPushButton(self)
        self.button_exit.setText("Exit")
        self.button_exit.setGeometry(x-60-60, y-40+m, 60, 30)
        self.button_exit.setStyleSheet("font-family: Arial;background-color: #e02424; color: white; font-size: 16px;")
        self.button_exit.clicked.connect(self.exit_program)
        main_layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout2 = QHBoxLayout()
        layout3 = QHBoxLayout()

        layout1.addWidget(self.label)
        layout1.addWidget(self.textbox)
        layout1.addWidget(self.button)

        layout2.addWidget(self.out_textbox)

        layout3.addWidget(self.button_clean)
        layout3.addWidget(self.button_exit)

        main_layout.addLayout(layout1)
        main_layout.addLayout(layout2)
        main_layout.addLayout(layout3)
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)
        # 下拉菜单栏目
        self.menubar = QMenuBar(self)

        tools = self.menubar.addMenu(QIcon("./imgs/tools.png"),"Tools")
        Onlinetools = self.menubar.addMenu("WebTools")
        helps = self.menubar.addMenu("Help")

        self.setStyleSheet("""
            QMenuBar {
                background-color: #666666;
                color: white;
            }
            QMenuBar::item {
                background-color: #666666;
                color: white;
                padding: 5px 10px;
                border: none;
            }
            QMenuBar::item:selected {
                background-color: #333;
            }

        """)
        self.setMenuBar(self.menubar)
        self.menubar.setFixedHeight(30)
        self.mass_title = "Mass fraction"
        self.mol_title = "Mol fraction"
        self.mdens_title = "Mass density"
        self.about_title = "About"
        self.online_molcalc_title = "MolCalc"
        self.online_moldraw_title = "MolDraw"
        self.online_molview_title = "MolView"
        massf = QAction(QIcon("./imgs/mass_men.png"),self.mass_title, self)
        molf  = QAction(QIcon("./imgs/mol_men.png"),self.mol_title, self)
        massd = QAction(QIcon("./imgs/dens_men.png"),self.mdens_title, self)
        about = QAction(QIcon("./imgs/about.png"),self.about_title, self)
        draw_mol = QAction(QIcon("./imgs/draw_mol.png"),self.online_moldraw_title, self)
        online_molcalc = QAction(QIcon("./imgs/online_molcalc.png"),self.online_molcalc_title, self)
        online_molview = QAction(QIcon("./imgs/online_molview.png"),self.online_molview_title, self)
        tools.addAction(massf)
        tools.addAction(molf)
        tools.addAction(massd)
        Onlinetools.addAction(draw_mol)
        Onlinetools.addAction(online_molview)
        Onlinetools.addAction(online_molcalc)

        helps.addAction(about)

        massf.triggered.connect(self.open_massfrac)
        molf.triggered.connect(self.open_molfrac)
        massd.triggered.connect(self.open_massdens)
        
        draw_mol.triggered.connect(self.open_DrawMol)
        online_molcalc.triggered.connect(self.open_MolCalc)
        online_molview.triggered.connect(self.open_MolView)

        about.triggered.connect(self.openAboutDialog)

    # about
    def openAboutDialog(self):
        dialog = AboutDialog()
        dialog.exec_()

    # draw molecules online
    def open_DrawMol(self):
        self.DrawMol = QWidget()
        self.DrawMol.setWindowTitle(self.online_molcalc_title)
        self.DrawMol.setWindowIcon(QIcon("./imgs/draw_mol.png"))
        self.DrawMol.setGeometry(200, 50, 1200, 900)
        self.DrawMol.setAttribute(Qt.WA_DeleteOnClose, False)

        layout = QVBoxLayout()
        widget = QWidget(self.DrawMol)
        widget.setLayout(layout)
        webview = QWebEngineView()
        layout.addWidget(webview)
        self.DrawMol.setLayout(layout)
        self.DrawMol.show()
        url = QUrl.fromUserInput("https://chemoinfo.ipmc.cnrs.fr/LEA3D/drawonline.html")  # 替换为你想要嵌入的网页的URL
        webview.load(url)
        webview.setZoomFactor(1.8)

    # online_molcalc
    def open_MolCalc(self):
        # self.molcalc = self.open_window("./imgs/online_molcalc.png",self.online_molcalc_title)
        self.molcalc = QWidget()
        self.molcalc.setWindowTitle(self.online_molcalc_title)
        self.molcalc.setWindowIcon(QIcon("./imgs/online_molcalc.png"))
        self.molcalc.setGeometry(200, 50, 1200, 900)
        self.molcalc.setAttribute(Qt.WA_DeleteOnClose, False)

        layout = QVBoxLayout()
        widget = QWidget(self.molcalc)
        widget.setLayout(layout)
        webview = QWebEngineView()
        layout.addWidget(webview)
        self.molcalc.setLayout(layout)
        self.molcalc.show()
        url = QUrl.fromUserInput("https://molcalc.org")  # 替换为你想要嵌入的网页的URL
        webview.load(url)

    # online_molview
    def open_MolView(self):
        self.molview = QWidget()
        self.molview.setWindowTitle(self.online_molview_title)
        self.molview.setWindowIcon(QIcon("./imgs/online_molview.png"))
        self.molview.setGeometry(200, 50, 1300, 900)
        self.molview.setAttribute(Qt.WA_DeleteOnClose, False)

        layout = QVBoxLayout()
        widget = QWidget(self.molview)
        widget.setLayout(layout)
        webview = QWebEngineView()
        layout.addWidget(webview)
        self.molview.setLayout(layout)
        self.molview.show()
        url = QUrl.fromUserInput("https://molview.org/")  # 替换为你想要嵌入的网页的URL
        webview.load(url)

    # mass fractions
    def open_massfrac(self):
        self.massfrac = self.open_window("./imgs/mass_men.png",self.mass_title)
        label = QtWidgets.QLabel(self.massfrac)
        label.setText("Formula:")
        label.move(20, 25)
        label.setStyleSheet("font-size: 18px;font-family: Arial;")   
        self.mass_textbox = QtWidgets.QLineEdit(self.massfrac)
        self.mass_textbox.setText("H2O 360 NaCl 4")
        self.mass_textbox.setGeometry(120, 20, 200, 30)
        self.mass_textbox.setStyleSheet("font-family: Arial;background-color: white; font-size: 18px;border: 1px solid #ccc;")

        button = QtWidgets.QPushButton(self.massfrac)
        button.setText("Run")
        button.setGeometry(350, 20, 55, 30)
        button.setStyleSheet("font-family: Arial;background-color: #808080; color: white; font-size: 18px;")
        button.clicked.connect(self.cal_massfrac)

        self.mass_out_textbox = QTextEdit(self.massfrac)
        self.mass_out_textbox.setGeometry(20, 70, 385, 140)
        self.mass_out_textbox.setReadOnly(True)
        self.mass_out_textbox.setStyleSheet("""
            QTextEdit {
                background-color: #f5f5f5;  /* 设置背景颜色 */
                font-family: Arial;  /* 设置字体样式 */
                font-size: 18px;
                border: 1px solid #ccc;  /* 设置边框 */
            }
        """)
        button_clean = QtWidgets.QPushButton(self.massfrac)
        button_clean.setText("Clean")
        button_clean.setGeometry(60, y-135, 60, 30)
        button_clean.setStyleSheet("font-family: Arial;background-color: #495c69; color: white; font-size: 18px;")
        button_clean.clicked.connect(self.clear_output_massfrac)

        button_exit = QtWidgets.QPushButton(self.massfrac)
        button_exit.setText("Exit")
        button_exit.setGeometry(x-60-60, y-135, 60, 30)
        button_exit.setStyleSheet("font-family: Arial;background-color: #e02424; color: white; font-size: 18px;")
        button_exit.clicked.connect(self.exit_massfrac)

        # ----------------------------------------------
        main_layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout2 = QHBoxLayout()
        layout3 = QHBoxLayout()

        layout1.addWidget(label)
        layout1.addWidget(self.mass_textbox)
        layout1.addWidget(button)

        layout2.addWidget(self.mass_out_textbox)

        layout3.addWidget(button_clean)
        layout3.addWidget(button_exit)

        main_layout.addLayout(layout1)
        main_layout.addLayout(layout2)
        main_layout.addLayout(layout3)
        self.massfrac.setLayout(main_layout)
        self.massfrac.show()
        # ----------------------------------------------

    def cal_massfrac(self):
        mol_formula_num = self.mass_textbox.text()
        self.mass_out_textbox.append(mol_formula_num)
        fraction_list = list(mol_formula_num.strip().split())
        try:
            molecular_mass_message_list = cf.mass_fraction(fraction_list)
            for message in molecular_mass_message_list:
                print(message)
                self.mass_out_textbox.append(str(message))
        except:
            message = "Input format error !!!"
            self.mass_out_textbox.append(message)
        return


    def open_molfrac(self):
        self.molfrac = self.open_window("./imgs/mol_men.png",self.mol_title)
        label = QtWidgets.QLabel(self.molfrac)
        label.setText("Formula:")
        label.move(20, 25)
        label.setStyleSheet("font-family: Arial;font-size: 18px;")   
        self.mol_textbox = QtWidgets.QLineEdit(self.molfrac)
        self.mol_textbox.setText("H2O 368 CH4 64")
        self.mol_textbox.setGeometry(120, 20, 200, 30)
        self.mol_textbox.setStyleSheet("font-family: Arial;background-color: white; font-size: 18px;border: 1px solid #ccc;")

        button = QtWidgets.QPushButton(self.molfrac)
        button.setText("Run")
        button.setGeometry(350, 20, 55, 30)
        button.setStyleSheet("font-family: Arial;background-color: #808080; color: white; font-size: 18px;")
        button.clicked.connect(self.cal_molfrac)

        self.mol_out_textbox = QTextEdit(self.molfrac)
        self.mol_out_textbox.setGeometry(20, 70, 385, 140)
        self.mol_out_textbox.setReadOnly(True)
        self.mol_out_textbox.setStyleSheet("""
            QTextEdit {
                background-color: #f5f5f5;  /* 设置背景颜色 */
                font-family: Arial;  /* 设置字体样式 */
                font-size: 18px;
                border: 1px solid #ccc;  /* 设置边框 */
            }
        """)
        button_clean = QtWidgets.QPushButton(self.molfrac)
        button_clean.setText("Clean")
        button_clean.setGeometry(60, y-135, 60, 30)
        button_clean.setStyleSheet("font-family: Arial;background-color: #495c69; color: white; font-size: 18px;")
        button_clean.clicked.connect(self.clear_output_molfrac)

        button_exit = QtWidgets.QPushButton(self.molfrac)
        button_exit.setText("Exit")
        button_exit.setGeometry(x-60-60, y-135, 60, 30)
        button_exit.setStyleSheet("font-family: Arial;background-color: #e02424; color: white; font-size: 18px;")
        button_exit.clicked.connect(self.exit_molfrac)
   
        # ----------------------------------------------
        main_layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout2 = QHBoxLayout()
        layout3 = QHBoxLayout()

        layout1.addWidget(label)
        layout1.addWidget(self.mol_textbox)
        layout1.addWidget(button)

        layout2.addWidget(self.mol_out_textbox)

        layout3.addWidget(button_clean)
        layout3.addWidget(button_exit)

        main_layout.addLayout(layout1)
        main_layout.addLayout(layout2)
        main_layout.addLayout(layout3)

        self.molfrac.setLayout(main_layout)

        self.molfrac.show()
        # ----------------------------------------------


    def cal_molfrac(self):
        mol_formula_num = self.mol_textbox.text()
        self.mol_out_textbox.append(mol_formula_num)
        fraction_list = list(mol_formula_num.strip().split())
        try:
            molecular_mol_message_list = cf.mole_fraction(fraction_list)

            for message in molecular_mol_message_list:
                self.mol_out_textbox.append(message)
        except:
            message = "Input format error!!!!"
            self.mol_out_textbox.append(message)

        return



    def open_massdens(self):
        self.massdens = self.open_window("./imgs/dens_men.png",self.mdens_title)
        label = QtWidgets.QLabel(self.massdens)
        label.setText("Formula:")
        label.move(20, 25)
        label.setStyleSheet("font-size: 18px;font-family: Arial;")   
        self.mass_dens_textbox = QtWidgets.QLineEdit(self.massdens)
        self.mass_dens_textbox.setText("H2O 7200")
        self.mass_dens_textbox.setGeometry(120, 20, 200, 30)
        self.mass_dens_textbox.setStyleSheet("font-family: Arial;background-color: white; font-size: 18px;border: 1px solid #ccc;")

        label_box = QtWidgets.QLabel(self.massdens)
        label_box.setText("Box:")
        label_box.setStyleSheet("font-size: 18px;font-family: Arial;")   

        self.mass_dens_textbox_box = QtWidgets.QLineEdit(self.massdens)
        self.mass_dens_textbox_box.setText("60 60 60")
        self.mass_dens_textbox_box.setGeometry(120, 20, 200, 30)
        self.mass_dens_textbox_box.setStyleSheet("font-family: Arial;background-color: white; font-size: 18px;border: 1px solid #ccc;")


        button = QtWidgets.QPushButton(self.massdens)
        button.setText("Run")
        button.setGeometry(350, 20, 55, 30)
        button.setStyleSheet("font-family: Arial;background-color: #808080; color: white; font-size: 18px;")
        button.clicked.connect(self.cal_massdens)

        self.mass_dens_out_textbox = QTextEdit(self.massdens)
        self.mass_dens_out_textbox.setGeometry(20, 70, 385, 140)
        self.mass_dens_out_textbox.setReadOnly(True)
        self.mass_dens_out_textbox.setStyleSheet("""
            QTextEdit {
                background-color: #f5f5f5;  /* 设置背景颜色 */
                font-family: Arial;  /* 设置字体样式 */
                font-size: 18px;
                border: 1px solid #ccc;  /* 设置边框 */
            }
        """)
        button_clean = QtWidgets.QPushButton(self.massdens)
        button_clean.setText("Clean")
        button_clean.setGeometry(60, y-135, 60, 30)
        button_clean.setStyleSheet("font-family: Arial;background-color: #495c69; color: white; font-size: 18px;")
        button_clean.clicked.connect(self.clear_output_massdens)

        button_exit = QtWidgets.QPushButton(self.massdens)
        button_exit.setText("Exit")
        button_exit.setGeometry(x-60-60, y-135, 60, 30)
        button_exit.setStyleSheet("font-family: Arial;background-color: #e02424; color: white; font-size: 18px;")
        button_exit.clicked.connect(self.exit_massdens)

        # ----------------------------------------------
        main_layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout2 = QHBoxLayout()
        layout3 = QHBoxLayout()

        layout1.addWidget(label)
        layout1.addWidget(self.mass_dens_textbox)
        layout1.addWidget(label_box)
        layout1.addWidget(self.mass_dens_textbox_box)
        layout1.addWidget(button)

        layout2.addWidget(self.mass_dens_out_textbox)

        layout3.addWidget(button_clean)
        layout3.addWidget(button_exit)

        main_layout.addLayout(layout1)
        main_layout.addLayout(layout2)
        main_layout.addLayout(layout3)
        self.massdens.setLayout(main_layout)
        self.massdens.show()
        # ----------------------------------------------
    def cal_massdens(self):
        mol_formula_num = self.mass_dens_textbox.text()
        xyz = self.mass_dens_textbox_box.text()
        self.mass_dens_out_textbox.append("Box Size = ("+xyz+") Å")
        self.mass_dens_out_textbox.append(mol_formula_num)
        fraction_list = list(mol_formula_num.strip().split())
        try:
            x,y,z = map(float,list(xyz.strip().split()))
            molecular_mass_message_list = cf.mass_density(fraction_list,x,y,z)
            self.mass_dens_out_textbox.append(molecular_mass_message_list)
            print(x,y,z)
        except:
            self.mass_dens_out_textbox.append("Input error!!!")

        return


    def clear_output_textbox(self):
    	self.out_textbox.clear()
    def clear_output_molfrac(self):
        self.mol_out_textbox.clear()
    def clear_output_massfrac(self):
        self.mass_out_textbox.clear()
    def clear_output_massdens(self):
        self.mass_dens_out_textbox.clear()

    def exit_program(self):
    	self.close()
    def exit_molfrac(self):
        self.molfrac.close()
    def exit_massfrac(self):
        self.massfrac.close()
    def exit_massdens(self):
        self.massdens.close()

    def calculate_molecular_weight(self):
        formula = self.textbox.text()
        m=cf.MoleculeMass()
        molmass = m.MolMass(formula)
        molmass = round(molmass,6)
        self.out_textbox.append("Formula = " + formula)
        self.out_textbox.append("Molecular Mass = " + str(molmass)+ " (g/mol)")

    def open_window(self,icon,title):
        window = QWidget()
        window.setWindowTitle(title)
        window.setWindowIcon(QIcon(icon))
        window.setGeometry(200, 200, x, y-100)
        # window.setMinimumSize(x, y-100)
        # window.setMaximumSize(mx, my)
        window.setAttribute(Qt.WA_DeleteOnClose, False)

        return window



if __name__ == "__main__":

    app = QApplication(sys.argv)
    window = ChemicalCalculator()
    window.show()
    sys.exit(app.exec_())                      