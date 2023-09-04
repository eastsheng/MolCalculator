#!/usr/bin/python
"""
calculating mole fraction X
n = m/M
n_fraction = m_fraction/M_fraction
n_total = m_total/M_total
X = n_fraction/n_total
n_fraction = Number of fractions
n_total = total Number

"""
import argparse
import sys
import re
import periodictable as pt

class MoleculeMass(object):
    """docstring for MassVolDen"""
    def __init__(self, ):
        super(MoleculeMass, self).__init__()
        # self.arg = arg

    def unnull(self,num):
        if num==[]:
            num=1
        else:
            num = int(num[0])
        return num  

    def Chem2Element(self,formula):
        '''chemical formula to element and number'''
        # print(formula)
        element = []
        number = []
        formula_sub = re.findall(r'[A-Z][^A-Z]*',formula)
        # print(formula_sub)

        for i in range(len(formula_sub)):
            # print(formula_sub[i])
            num = re.findall(r'\d+',formula_sub[i])
            num = self.unnull(num)
            # print(num)
            number.append(num)
            ele = ''.join(re.findall(r'[A-Za-z]',formula_sub[i]))
            ele = ele.title()
            element.append(ele)
        # print(element,number)

        return element, number

    def MolMass(self,formula):
        '''分子质量，给出分子中各元素的个数，返回单个分子的质量（g）'''
        element, number = self.Chem2Element(formula)
        molmass = 0
        
        for i in range(len(element)):
            try:
                ele_mass = pt.elements.symbol(element[i]).mass
            except:
                print("ERROR: Can't find "+element[i]+" ......checkout......")
                break     
            molmass += ele_mass*number[i]
        # print(formula,"Molecular mass = ",molmass,"g/mol")
        
        return molmass


def calculate_molar_fraction(compound_dict):
    total_amount = sum(compound_dict.values())
    molar_fractions = {}

    for compound, amount in compound_dict.items():
        molar_fraction = amount / total_amount
        molar_fractions[compound] = molar_fraction

    return molar_fractions

def list_to_dict(lst):
    # 使用字典推导式将列表转换为字典
    # 偶数索引位置的元素作为键，奇数索引位置的元素作为对应键的值
    dictionary = {lst[i]: int(lst[i+1]) for i in range(0, len(lst), 2)}
    return dictionary


def mole_fraction(molecules_list):
    compound_dict = list_to_dict(molecules_list)
    # 传入化合物及其数量的字典作为参数
    # compound_dict = {"H2O": 368, "CH4": 64}

    # 调用函数计算摩尔分数
    result = calculate_molar_fraction(compound_dict)

    # 打印结果
    result_list = []
    print("\n#",20*"-","Calculating mole fraction",20*"-","#\n")
    for compound, molar_fraction in result.items():
        print("Mole fraction of",compound, "=", round(molar_fraction,6))
        result_list.append("Mole fraction of "+compound+" = "+str(round(molar_fraction,6)))
    print("\n#",20*"-","Mole fraction end!!!!!!!!",20*"-","#\n")
    

    return result_list    


def cal_total_mass(molecules_list):
    mol_dict = list_to_dict(molecules_list)
    # print(mol_dict)
    m=MoleculeMass()
    mol_symbol = list(mol_dict.keys())
    mol_number = list(mol_dict.values())
    total_mass = 0
    for i in range(len(mol_dict)):
        molecular_mass = m.MolMass(mol_symbol[i]) 
        total_mass += molecular_mass*mol_number[i]
    total_mass_message = "Total mass of molecules ="+str(total_mass)    
    print(total_mass_message)
    return total_mass 

def mass_fraction(molecules_list):
    print("\n#",20*"-","Claculating mass fraction",20*"-","#\n")
    total_mass_message = cal_total_mass(molecules_list)
    molecular_mass_message_list = []
    molecular_mass_message_list.append(total_mass_message)
    for i in range(len(mol_dict)):
        molecular_mass = m.MolMass(mol_symbol[i]) 
        mf = molecular_mass*mol_number[i]/total_mass
        molecular_mass_message = "Mass fraction of "+str(mol_symbol[i])+" = "+str(round(mf,6))
        # return molecular_mass_message
        print(molecular_mass_message)
        molecular_mass_message_list.append(molecular_mass_message)
    # print(molecular_mass_message_list)
    print("\n#",20*"-","Mass fraction end!!!!!!!!",20*"-","#\n")
    # molecular_mass_message_list.append(mass_end_message)
    return molecular_mass_message_list


amu2g = 6.02214076208112e23
A2cm = 1e-8

def mass_density(molecules_list,x=22,y=22,z=22):
    '''给出总质量和三个方向尺寸(埃米)，返回密度（g/mL）'''
    tot_mass = cal_total_mass(molecules_list)
    unitconvert = amu2g*(A2cm)**3
    vol = (x*y*z)*unitconvert
    dens = tot_mass/vol
    dens = round(dens,6)
    dens_out = "Mass Density = "+str(dens)+" (g/mL)"
    # print(dens_out)
    return dens_out



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fraction calculation functions')
    parser.add_argument('-mole', '--molefrac', type=str, 
        help='name and number of fractions,for example: -mole "H2O 368 CH4 64"')
    parser.add_argument('-mass', '--massfrac', type=str, 
        help='formula and number of fractions,for example: -mass "H2O 360 NaCl 4"')
    
    args = parser.parse_args()
    print('\nCommand line:')
    print("\t\t$python"," ".join(sys.argv))

    """calculating mole fraction"""
    try:
        fraction_list = list(args.molefrac.strip().split())
        mole_fraction(fraction_list)
    except:
        pass

    """calculating mass fraction"""
    # # mol_formula_num = "H2O 368 CH4 64"
    # mol_formula_num = args.massfrac
    try:
        mol_list = list(args.massfrac.strip().split())
        mass_fraction(mol_list)
    except:
        pass