################################################################
#                                                              #
# LOVCI: Ladder Operator Vibrational Configuration Interaction #
#                                                              #
################################################################

#LOVCI semi-automated test suite

#Modules
import subprocess
import sys
import os

#Classes
class ClrSet:
  #Unicode colors
  Norm = '\033[0m'
  Red = '\033[91m'
  Bold = '\033[1m'
  Green = '\033[92m'
  Blue = '\033[94m'
  Cyan = '\033[36m'
  #Set colors
  TFail = Bold+Red #Highlight failed tests
  TPass = Bold+Green #Highlight passed tests
  Reset = Norm #Reset to defaults

#Print title
line = '\n'
line += "****************************************************************"
line += '\n'
line += "*                                                              *"
line += '\n'
line += "* LOVCI: Ladder Operator Vibrational Configuration Interaction *"
line += '\n'
line += "*                                                              *"
line += '\n'
line += "****************************************************************"
line += '\n'
print(line)

#Check for help
if (len(sys.argv) >= 2):
  #Print help if arguments are given
  line = ""
  line += "Usage:"
  line += '\n'
  line += " user:$ ./runtests"
  line += '\n'
  print(line)

#Extract CPUs
try:
  cmd = "echo ${OMP_NUM_THREADS}"
  Ncpus = subprocess.check_output(cmd,shell=True)
  Ncpus = Ncpus.strip() #Saved as a string
except:
  Ncpus = "1"

#Print output
result = "Running tests with "
result += Ncpus
result += " CPUs..."
result += '\n'
result += '\n'
result += "Test results:"
print(result)

#First test (harmonic)
os.chdir("Harmonic")
result = ""
ZPE = ""
RunTime = ""
cmd = "lovci -n "+Ncpus+" -i Harm.inp -o Spect.txt > Log.txt"
subprocess.call(cmd,shell=True)
ifile = open("Log.txt")
data = ifile.readlines()
ifile.close()
for line in data:
  line = line.strip()
  line = line.split()
  if (len(line) > 3):
    if (line[0] == "Zero-point"):
      ZPE = line[2]
    if ((line[0] == "Run") and (line[1] == "time:")):
      RunTime = line[2]+" "+line[3]
data = []
result += " Harmonic oscillator: "
if (ZPE == "2825.00"):
  result += ClrSet.TPass+"Pass"+ClrSet.Reset+", "
else:
  result += ClrSet.TFail+"Fail"+ClrSet.Reset+", "
result += RunTime
print(result)
cmd = "rm -f Spect.txt Log.txt"
subprocess.call(cmd,shell=True)
os.chdir("../")

#Second test (1D cubic)
os.chdir("1D_cubic")
result = ""
ZPE = ""
RunTime = ""
cmd = "lovci -n "+Ncpus+" -i Cubic.inp -o Spect.txt > Log.txt"
subprocess.call(cmd,shell=True)
ifile = open("Log.txt")
data = ifile.readlines()
ifile.close()
for line in data:
  line = line.strip()
  line = line.split()
  if (len(line) > 3):
    if (line[0] == "Zero-point"):
      ZPE = line[2]
    if ((line[0] == "Run") and (line[1] == "time:")):
      RunTime = line[2]+" "+line[3]
data = []
result += " Anharmonic oscillator: "
if (ZPE == "497.57"):
  result += ClrSet.TPass+"Pass"+ClrSet.Reset+", "
else: 
  result += ClrSet.TFail+"Fail"+ClrSet.Reset+", "
result += RunTime
print(result)
cmd = "rm -f Spect.txt Log.txt"
subprocess.call(cmd,shell=True)
os.chdir("../")

#Third test (C2v potential)
os.chdir("C2v_pot")
result = ""
ZPE = ""
RunTime = ""
cmd = "lovci -n "+Ncpus+" -i C2v.inp -o Spect.txt > Log.txt"
subprocess.call(cmd,shell=True)
ifile = open("Log.txt")
data = ifile.readlines()
ifile.close()
for line in data:
  line = line.strip()
  line = line.split()
  if (len(line) > 3):
    if (line[0] == "Zero-point"):
      ZPE = line[2]
    if ((line[0] == "Run") and (line[1] == "time:")):
      RunTime = line[2]+" "+line[3]
data = []
result += " Four mode progression: "
if (ZPE == "4374.59"):
  result += ClrSet.TPass+"Pass"+ClrSet.Reset+", "
else:
  result += ClrSet.TFail+"Fail"+ClrSet.Reset+", "
result += RunTime
print(result)
cmd = "rm -f Spect.txt Log.txt"
subprocess.call(cmd,shell=True)
os.chdir("../")

#Print results
result = '\n'
result += "Done."
result += '\n'
print(result)
