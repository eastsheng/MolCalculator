"""
# 打包
pyinstaller -D -w ./MolCalculator.py --clean -i ./imgs/icons-64.png
# 一个分子计算器，可以实现：
# 1. 输入化学式，输出摩尔质量
# 2. 输入分子化学式和分子数目（至少两组），输出质量分数
# 3. 输入分子化学式和分子数目（至少两组），输出摩尔分数
# 4. 输入分子化学式和分子数目，与体系大小(Å)，输出质量密度g/mL
# 5. 添加Draw molecules, MolCalc和Molview三个WEB工具
# 6. 添加元素周期表WEB网页
# 7. 添加通过名称或者SMILES码搜索分子信息
"""
import sys
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QHBoxLayout ,QVBoxLayout,QTabWidget,\
QTextEdit, QMenuBar, QMenu, QAction, QMessageBox, QWidget, QVBoxLayout, QLabel, QDialog, QComboBox
from PyQt5.QtGui import QColor, QFont, QIcon, QPixmap
from PyQt5.QtCore import Qt, QCoreApplication, QSize, QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView
from qt_material import apply_stylesheet

import scripts.cal_fraction as cf
import scripts.FindCompounds as fc

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw


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

        # tools = self.menubar.addMenu(QIcon("./imgs/tools.png"),"Tools")
        tools = self.menubar.addMenu("🛠Tools")
        Onlinetools = self.menubar.addMenu("🌍WebTools")
        SearchMols = self.menubar.addMenu("🔍SearchMols")
        helps = self.menubar.addMenu("💡Help")

        self.setMenuBar(self.menubar)
        self.menubar.setFixedHeight(30)
        self.mass_title = "Mass fraction"
        self.mol_title = "Mol fraction"
        self.mdens_title = "Mass density"
        self.setting_title = "Setting"
        self.about_title = "About"
        self.PeriodicTable_title = "PeriodicTable"
        self.online_molcalc_title = "MolCalc"
        self.online_moldraw_title = "MolDraw"
        self.online_molview_title = "MolView"
        self.findmols_title = "FindMols"
        massf = QAction(QIcon("./imgs/mass_men.png"),self.mass_title, self)
        molf  = QAction(QIcon("./imgs/mol_men.png"),self.mol_title, self)
        massd = QAction(QIcon("./imgs/dens_men.png"),self.mdens_title, self)
        setting = QAction(QIcon("./imgs/setting.png"),self.setting_title, self)
        about = QAction(QIcon("./imgs/about.png"),self.about_title, self)
        periodictable = QAction(QIcon("./imgs/periodictable.png"),self.PeriodicTable_title, self)
        draw_mol = QAction(QIcon("./imgs/draw_mol.png"),self.online_moldraw_title, self)
        online_molcalc = QAction(QIcon("./imgs/online_molcalc.png"),self.online_molcalc_title, self)
        online_molview = QAction(QIcon("./imgs/online_molview.png"),self.online_molview_title, self)
        
        findmols = QAction(QIcon("./imgs/find_mols.png"),self.findmols_title, self)
        self.setStyleSheet("""
            QMenuBar {
                background-color: #666666;
                color: white;
                font-family: Arial;
                font: bold;
            }
            QMenuBar::item {
                background-color: #666666;
                color: white;
                font-size: 14px;
                padding: 5.5px 10px;

            }
            QMenuBar::item:selected {
                background-color: #333;
            }
        """)        

        tools.addAction(massf)
        tools.addAction(molf)
        tools.addAction(massd)
        Onlinetools.addAction(periodictable)
        Onlinetools.addAction(draw_mol)
        Onlinetools.addAction(online_molview)
        Onlinetools.addAction(online_molcalc)
        SearchMols.addAction(findmols)

        helps.addAction(setting)
        helps.addAction(about)

        massf.triggered.connect(self.open_massfrac)
        molf.triggered.connect(self.open_molfrac)
        massd.triggered.connect(self.open_massdens)
        
        periodictable.triggered.connect(self.open_periodictable)
        draw_mol.triggered.connect(self.open_DrawMol)
        online_molcalc.triggered.connect(self.open_MolCalc)
        online_molview.triggered.connect(self.open_MolView)

        findmols.triggered.connect(self.open_FindMols)

        setting.triggered.connect(self.openSetting)
        about.triggered.connect(self.openAboutDialog)

    # about
    def openAboutDialog(self):
        dialog = AboutDialog()
        dialog.exec_()

    # setting
    def openSetting(self):
        self.Setting = self.open_window("./imgs/setting.png",self.setting_title,y=100)
        # self.Setting.setMaximumSize(mx, 30)

        layout = QVBoxLayout(self.Setting)
        tab_widget = QTabWidget()
        layout.addWidget(tab_widget)

        themes = QWidget()
        tab_widget.addTab(themes, "🎨 Themes")
        themes_layout = QHBoxLayout(themes)
        themes_layout.addWidget(QLabel("Select theme (Restart)："))
        tab_widget.setStyleSheet("font-family: Arial; font-size: 15px; font: bold;")

        themes = QWidget()
        tab_widget.addTab(themes, "💤 Others")

        
        self.theme_combo = QComboBox()
        self.themes_names_dict = read_themes()
        self.themes_names = list(self.themes_names_dict.keys())
        self.theme_combo.addItems(self.themes_names)
        self.theme_combo.setCurrentText(default_theme)
        themes_layout.addWidget(self.theme_combo)

        self.theme_combo.currentIndexChanged[str].connect(self.write_chosen_theme)
        
        self.theme_combo.setStyleSheet(
            "font-family: Arial; font-size: 15px; font: bold; background-color: #f5f5f5;")
        self.Setting.setLayout(layout)
        self.Setting.show()

    def write_chosen_theme(self):
        try:
            # if self.theme_combo.currentText() in self.themes_names:
            choose_theme_name = self.theme_combo.currentText()
            choose_theme=self.themes_names_dict[choose_theme_name]
            with open(chosen_theme_file,"w") as t:
                t.write(choose_theme)
                print(choose_theme,"theme was selected......")
        except:
            with open(chosen_theme_file,"w") as t:
                t.write(default_theme) 
                print("default "+default_theme+" theme was selected......")

        return


    # open_FindMols
    def open_FindMols(self):
        self.FindMols = self.open_window("./imgs/find_mols.png",self.findmols_title,y=600)
        label = QtWidgets.QLabel(self.FindMols)
        label.setText("Name or SMILES:")
        label.move(20, 25)
        label.setStyleSheet("font-size: 18px;font-family: Arial;")  
        self.display_label = QtWidgets.QLabel(self.FindMols) 
        # self.display_label.setText("Molecular Picture:")
        self.display_label.setFixedSize(x,200)
        self.display_label.setStyleSheet("font-size: 18px;font-family: Arial;background-color: white;")  
        self.display_label.setAlignment(Qt.AlignCenter)
        self.findmols_textbox = QtWidgets.QLineEdit(self.FindMols)
        self.findmols_textbox.setText("PVP")
        self.findmols_textbox.setGeometry(120, 20, 200, 30)
        self.findmols_textbox.setStyleSheet("font-family: Arial;background-color: white; font-size: 18px;border: 1px solid #ccc;")

        button = QtWidgets.QPushButton(self.FindMols)
        button.setText("Run")
        button.setGeometry(350, 20, 55, 30)
        button.setStyleSheet("font-family: Arial;background-color: #808080; color: white; font-size: 18px;")
        
        button.clicked.connect(self.getMols)

        self.findmols_out_textbox = QTextEdit(self.FindMols)
        self.findmols_out_textbox.setGeometry(20, 70, 385, 140)
        self.findmols_out_textbox.setReadOnly(True)
        self.findmols_out_textbox.setStyleSheet("""
            QTextEdit {
                background-color: white;  /* 设置背景颜色 */
                font-family: Arial;  /* 设置字体样式 */
                font-size: 18px;
                border: 0.1px solid #ccc;  /* 设置边框 */
            }
        """)
        button_clean = QtWidgets.QPushButton(self.FindMols)
        button_clean.setText("Clean")
        button_clean.setGeometry(60, y-135, 60, 30)
        button_clean.setStyleSheet("font-family: Arial;background-color: #495c69; color: white; font-size: 18px;")
        button_clean.clicked.connect(self.clear_output_findmols)

        button_exit = QtWidgets.QPushButton(self.FindMols)
        button_exit.setText("Exit")
        button_exit.setGeometry(x-60-60, y-135, 60, 30)
        button_exit.setStyleSheet("font-family: Arial;background-color: #e02424; color: white; font-size: 18px;")
        button_exit.clicked.connect(self.exit_findmols)

        # ----------------------------------------------
        main_layout = QVBoxLayout()
        layout1 = QHBoxLayout()
        layout12 = QHBoxLayout()
        layout2 = QHBoxLayout()
        layout3 = QHBoxLayout()

        layout1.addWidget(label)
        layout1.addWidget(self.findmols_textbox)
        layout1.addWidget(button)
        layout12.addWidget(self.display_label)
        layout2.addWidget(self.findmols_out_textbox)

        layout3.addWidget(button_clean)
        layout3.addWidget(button_exit)

        main_layout.addLayout(layout1)
        main_layout.addLayout(layout2)
        main_layout.addLayout(layout12)
        main_layout.addLayout(layout3)
        self.FindMols.setLayout(main_layout)
        self.FindMols.show()
        # ----------------------------------------------

    def getMols(self):
        name = self.findmols_textbox.text()
        compounds = fc.get_type_compounds(name=name)
        mols_list,imgs_list,infos_list = fc.mol_infos(compounds)
        n = len(mols_list)
        if n == 0:
            self.findmols_out_textbox.append("Warning: Not found '"
                +name+"',there is no 'Name' or 'SMILES' in 'Pubchem' database")
        else:       
            for i in range(n):
                imgs_list[i].save('./temp/temp.png')
                image = QPixmap('./temp/temp.png').scaled(200, 200)
                self.display_label.setPixmap(image)
                self.findmols_out_textbox.append("SMILES: "+str(infos_list[i]["smiles"]))
                self.findmols_out_textbox.append("Formula: "+str(infos_list[i]["formula"]))
                self.findmols_out_textbox.append("Weight: "+str(infos_list[i]["weight"]))
                self.findmols_out_textbox.append("IUPAC NAME: "+str(infos_list[i]["name"]))
                self.findmols_out_textbox.append("\n")


    # online Periodic Table
    def open_periodictable(self):
        self.periodictable = QWidget()
        self.periodictable.setWindowTitle(self.PeriodicTable_title)
        self.periodictable.setWindowIcon(QIcon("./imgs/periodictable.png"))
        self.periodictable.setGeometry(200, 50, 1200, 900)
        self.periodictable.setAttribute(Qt.WA_DeleteOnClose, False)

        layout = QVBoxLayout()
        widget = QWidget(self.periodictable)
        widget.setLayout(layout)
        webview = QWebEngineView()
        layout.addWidget(webview)
        self.periodictable.setLayout(layout)
        self.periodictable.show()
        url = QUrl.fromUserInput("https://www.rsc.org/periodic-table")  # 替换为你想要嵌入的网页的URL
        webview.load(url)


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

    def clear_output_findmols(self):
        self.findmols_out_textbox.clear()

    def clear_output_molfrac(self):
        self.mol_out_textbox.clear()
    def clear_output_massfrac(self):
        self.mass_out_textbox.clear()
    def clear_output_massdens(self):
        self.mass_dens_out_textbox.clear()

    def exit_program(self):
    	self.close()
    def exit_findmols(self):
        self.FindMols.close()
    
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

    def open_window(self,icon,title,y=y):
        window = QWidget()
        window.setWindowTitle(title)
        window.setWindowIcon(QIcon(icon))
        window.setGeometry(200, 200, x, y-100)
        # window.setMinimumSize(x, y-100)
        # window.setMaximumSize(mx, my)
        window.setAttribute(Qt.WA_DeleteOnClose, False)

        return window

def read_themes():
    with open(all_themes_file,"r") as t:
        themes = t.readlines()
    # print(themes)
    themes_names = {}
    names,files = [],[]
    for theme in themes:
        name = theme.split(".")[0]
        names.append(name) 
        files.append(theme.split("\n")[0]) 
    pairs = list(zip(names, files))
    themes_names = dict(pairs)
    return themes_names

def read_choose_theme():
    try:
        with open(chosen_theme_file,"r") as t:
            chosen_theme = t.readlines()[0].split("\n")[0]
            print(chosen_theme,"theme was selected......")
    except:
        print("Warning: Not found chosen theme......")
        print("Warning: default "+default_theme+" theme was selected......")
        chosen_theme = default_theme
    
    return chosen_theme

if __name__ == "__main__":
    all_themes_file = "./scripts/themes"
    chosen_theme_file = "./scripts/.theme"

    default_theme = "light_blue.xml"
    default_theme = default_theme.split(".")[0]
    app = QApplication(sys.argv)
    # setup stylesheet
    chosen_theme = read_choose_theme()
    # theme = themes_names["dark_blue"]
    apply_stylesheet(app, chosen_theme)
    window = ChemicalCalculator()
    window.show()
    sys.exit(app.exec_())                      