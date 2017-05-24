#coding=utf-8
import os, sys, re
import numpy as np


Bohr_to_An = 0.52917721092


def read_gmx_gjf(inp):
	#return data of coord and charge
	with open(inp, "r") as f:
		inp_text = f.readlines()

	coord, charge = [], []
	iscoord = True
	for line in inp_text:
		if not line.strip() and iscoord and len(coord):
			iscoord = False
		elif line.strip() and line.strip()[0].isnumeric() and len([i for i in line.split(' ') if i]) == 4:
			if iscoord:
				coord.append(line)
			else:
				charge.append(line)

	for num in range(len(coord)):
		natom, x, y, z = [_ for _ in coord[num].strip().split(' ') if _]
		coord[num] = [int(natom), float(x) * Bohr_to_An, float(y) * Bohr_to_An, float(z) * Bohr_to_An]

	for num in range(len(charge)):
		x, y, z, chrg = [_ for _ in charge[num].strip().split(' ') if _]
		charge[num] = [float(x) * Bohr_to_An, float(y) * Bohr_to_An, float(z) * Bohr_to_An, float(chrg)]

	return coord, charge

def write_gau_fort7(energy, force_list, field_list, charge):
	with open("fort.7", "w") as f:
		f.write("%f\n"%energy)
		for line in force_list:
			if len(line) == 1:
				continue
			x, y, z = line
			f.write("%20.10E %20.10E %20.10E\n"%(- x, - y, - z))
		for num, line in enumerate(field_list):
			x, y, z = line
			f.write("%20.10E %20.10E %20.10E\n"%(- x, - y, - z))
		f.write("\n")

class QuantumCalculateFail(BaseException):
	pass


class BaseCalculator(object):

	def __init__(self, sys_charge = 0, sys_mlt = 1, mid_name = 'inp.gjf', mid_return = "inp.out"):
		self.CHARGE = sys_charge
		self.MLT = sys_mlt
		self.MIDRUN = mid_name
		self.MIDRETURN = mid_return

	def _preprocessing(self):
		#delete tmp files
		pass

	def _execute(self):
		return ""

	def _mk_mid_file(self, coord, charge):
		template = ""
		return template

	def calculate(self, coord, charge):
		self._preprocessing()
		with open(self.MIDRUN, "w") as f:
			f.write(self._mk_mid_file(coord, charge))
		os.system(self._execute())
		return self._read_return()

	def _read_return(self):
		#Fulllfill energy, force_list and field_list
		energy = 0.0
		force_list, field_list = [], []
		self.energy, self.force_list, self.field_list = energy, force_list, field_list
		return energy, force_list, field_list


class GaussCalculator(BaseCalculator):
	"""
	
	"""

	def __init__(self, system = "Linux", gauss_exe = "g09", scfcyc = 128, 
		method = "HF", basis = "6-31G", nproc = 1, mem = "200MW", 
		sys_charge = 0, sys_mlt = 1, mid_name = 'inp.gjf', mid_return = "inp.out"):
		super(GaussCalculator, self).__init__(sys_charge = sys_charge, sys_mlt = sys_mlt, mid_name = mid_name, mid_return = mid_return)
		self.SYSTEM = system
		self.PROGRAME = gauss_exe
		self.NPROC = nproc
		self.MEM = mem
		self.COMMAND = "# {method}/{basis} FORCE Nosymm Punch=Derivatives".format(method = method, basis = basis)
		self.RETRY = 0

	def _execute(self):
		if self.SYSTEM == "Linux":
			return "{GAU} < {MIDRUN} > {MIDRETURN}".format(GAU=self.PROGRAME, MIDRUN=self.MIDRUN, MIDRETURN=self.MIDRETURN)
		else:
			return "{GAU} {MIDRUN}".format(GAU=self.PROGRAME, MIDRUN=self.MIDRUN)

	def _preprocessing(self, delCHK = False):
		listdir = os.listdir('.')
		if "fort.7" in listdir:
			os.remove("fort.7")
		for fname in listdir:
			if "gxx." in fname or "Gau" in fname:
				os.remove(fname)
			if delCHK and ".chk" in fname:
				os.remove(fname)

	def _mk_mid_file(self, coord, charge):
		self._inp_coord, self._inp_charge = coord, charge
		_COMMAND = self.COMMAND
		if len(charge):
			_COMMAND += " CHARGE Prop=(Field, Read)"
		if "inp.chk" in os.listdir('.'):
			_COMMAND += " GUESS=READ"

		template = ""
		template += "%MEM={MEM}\n%CHK=inp.chk\n%NPROC={NPROC}\n".format(MEM=self.MEM, NPROC=self.NPROC)

		template += _COMMAND + "\n\n"

		template += "GMX_INP\n\n{CHRG} {MLT}\n".format(CHRG=self.CHARGE, MLT=self.MLT)

		for natom, x, y, z in coord:
			template += "%i %16.10f %16.10f %16.10f \n"%(natom, x, y, z)
		template += "\n"
		for x, y, z, chrg in charge:
			template += "%16.10f %16.10f %16.10f %16.10f \n"%(x, y, z, chrg)
		template += "\n"
		for x, y, z, chrg in charge:
			template += "%16.10f %16.10f %16.10f \n"%(x, y, z)
		template += "\n\n"

		return template

	def _read_return(self):
		try:
			with open(self.MIDRETURN, "r") as f:
				out_text = f.readlines()
			with open("fort.7") as f:
				force_list = f.readlines()[1:]
		except FileNotFoundError as e:
			if self.RETRY > 2:
				raise e
			self.RETRY += 1
			self._preprocessing(delCHK=True)
			return self.calculate(self._inp_coord, self._inp_charge)

		energy = None
		for line in out_text:
			if "SCF Done" in line:
				energy = float(re.findall(r"[-0-9.]{9,}", line)[0])
		
		force_list = [[float('E'.join(j.split('D'))) for j in i.strip().split(' ') if j] for i in force_list if i.strip()]
		if not len(force_list):
			raise QuantumCalculateFail("No data in 'fort.7'.")

		for num in range(len(out_text)):
			if "-------- Electric Field --------" in out_text[num]:
				break
		start = num + 3 + len(self._inp_coord)
		end = num + 3 + len(self._inp_coord) + len(self._inp_charge)
		field_list = [[float(j) for j in i.strip().split(' ') if j][-3:] for i in out_text[start:end]]
		if len(field_list) != len(self._inp_charge):
			raise QuantumCalculateFail("Field length is only %i while charge length is %i"%(len(field_list), len(self._inp_charge)))

		self.energy, self.force_list, self.field_list = energy, force_list, field_list
		return energy, force_list, field_list


class DFTBPCalculator(BaseCalculator):

	atom_index = {1:['H', 's', '-0.1857'], 
			  6:['C', 'p', '-0.1492'], 
			  7:['N', 'p', '-0.1535'], 
			  8:['O', 'p', '-0.1575'], 
			  9:['F', 'p', '-0.1623'], 
			  15:['P', 'p', '-0.1400'], 
			  16:['S', 'd', '-0.1100'], 
			  17:['Cl', 'd', '-0.0697']}

	def __init__(self,prefix = None, num_iter = 128, sys_charge = 0, sys_mlt = 1, mid_name = 'dftb_in.hsd', mid_return = "detailed.out"):
		super(DFTBPCalculator, self).__init__(sys_charge = sys_charge, sys_mlt = sys_mlt, mid_name = mid_name, mid_return = mid_return)
		self.NUM_ITER = num_iter
		self.PREFIX = prefix
		if not self.PREFIX:
			raise QuantumCalculateFail("S-K parameters need to be set.")

	def _execute(self):
		return "dftb+ > output.log"

	def _preprocessing(self):
		listdir = os.listdir('.')
		for fname in listdir:
			if ".out" in fname or ".bin" in fname or ".hsd" in fname or "tmp" in fname or 'fort.7' in fname:
				os.remove(fname)

	def _mk_mid_file(self, coord, charge):
		self._inp_coord, self._inp_charge = coord, charge

		template = "Geometry = {\n"
		atom_to_num = {a:num for num, a in enumerate(set([i[0] for i in coord]))}
		num_to_atom = {v:k for k, v in atom_to_num.items()}
		template += """TypeNames = {"""
		for i in range(len(atom_to_num)):
			template += ' "' + self.atom_index[num_to_atom[i]][0] + '"'
		template += "}\n"
		template += "TypesAndCoordinates[Angstrom] = {\n"
		for natom, x, y, z in coord:
			template += "%i %16.8f %16.8f %16.8f\n"%(atom_to_num[natom]+1, x, y, z)
		template += "}\n}\n"
		template += "Hamiltonian = DFTB{\n"
		template += "SCC = Yes\nMaxSCCIterations = {num_iter}\nSCCTolerance = 1e-6\n".format(num_iter=self.NUM_ITER)
		template += "Filling = Fermi {\nTemperature [K] = 300\n}\n"
		template += "MaxAngularMomentum = {\n"
		for i in range(len(atom_to_num)):
			template += self.atom_index[num_to_atom[i]][0] + ' = "' + self.atom_index[num_to_atom[i]][1] + '"\n'
		template += "}\n"
		template += "Charge = {charge}\n".format(charge=self.CHARGE)
		template += "SlaterKosterFiles = Type2FileNames {\nPrefix = '%s'\nSeparator = '-'\nSuffix = '.skf'\n}\n"%self.PREFIX
		template += "Dispersion = DftD3{\nDamping = BeckeJohnson{\na1 = 0.746\na2 = 4.191\n}\ns8 = 3.209\n}\n"
		template += "ThirdOrderFull = Yes\nDampXH = Yes\nDampXHExponent = 4.00\nHubbardDerivs {\n"
		for i in range(len(atom_to_num)):
			template += self.atom_index[num_to_atom[i]][0] + ' = ' + self.atom_index[num_to_atom[i]][2] + '\n'
		template += '}\n'
		
		if len(charge):
			template += "ElectricField = {\n"
			template += "PointCharges = {\n"
			template += "CoordsAndCharges [Angstrom] = {\n"
			for x, y, z, chrg in charge:
				template += "%16.8f %16.8f %16.8f %f\n"%(x, y, z, chrg)
			template += "}\n}\n}\n"

		template += 'ForceEvaluation = "dynamics"\n}\nAnalysis = {\nCalculateForces = Yes\n}\n'

		return template

	def _read_return(self):
		with open(self.MIDRETURN, "r") as f:
			out_text = f.readlines()
		for line in out_text:
			if "Total energy" in line:
				energy = float([i for i in line.split(' ') if i][2]) #Hatree
		for num, line in enumerate(out_text):
			if "Total Forces" in line:
				break
		force_list = [[ float(j) for j in i.split(' ') if j] for i in out_text[num+1:num+1+len(self._inp_coord)]]

		if len(force_list) == 0:
			raise QuantumCalculateFail("SCC Error.")
		
		field_list = []
		if len(self._inp_charge):
			for num, line in enumerate(out_text):
				if "Forces on external charges" in line:
					break
			field_list = [[float(j) for j in i.split(' ') if j] for i in out_text[num+1:num+1+len(self._inp_charge)]]

			try:
				self.energy, self.force_list, self.field_list = energy, force_list, field_list
			except NameError as e:
				raise QuantumCalculateFail("Calculation Error.")
		return energy, force_list, field_list
