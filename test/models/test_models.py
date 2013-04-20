#!/usr/bin/env python

from subprocess import call
import os
import sys

ester_root=sys.argv[1]

def exec_cmd(cmd,infile=None,outfile=None):
	cmd=ester_root+"/bin/"+cmd
	if infile:
		f_in=open(infile,"r")
	else:
		f_in=None
	if outfile:
		f_out=open(outfile,"w")
	else:
		f_out=None
		
	try:
		out=call(cmd.split(),stdin=f_in,stdout=f_out)
	except:
		print("TEST FAILED !!!")
		sys.exit(1)
	if f_in:
		f_in.close()
	if f_out:
		f_out.close()
	if out:
		print("TEST FAILED !!!")
		sys.exit(1)
	return out

def compare_files(file1,file2):
	try:
		f1=open(file1,"r");
	except:
		return(False)
	try:
		f2=open(file2,"r");
	except:
		f1.close()
		return(False)
	if not f1.read()==f2.read():
		f1.close()
		f2.close()
		return False
	else:
		f1.close()
		f2.close()
		return True
			

print("Test model #1:")
cmd="ester 1d -M 5" 
cmd=cmd+" -p 1d.par -noplot -tol 1e-8 -maxit 100 -o test_model1"
exec_cmd(cmd)
cmd="ester output test_model1"
exec_cmd(cmd,"template_1d","test_out1")
test_result=compare_files("out1","test_out1")
if test_result:
	print("TEST OK\n")
else:
	print("TEST FAILED !!!")
	sys.exit(1)

print("Test model #2:")
cmd="ester 2d -i test_model1 -Omega_bk 0.5" 
cmd=cmd+" -p 2d.par -noplot -tol 1e-8 -maxit 10 -o test_model2"
exec_cmd(cmd)
cmd="ester output test_model2"
exec_cmd(cmd,"template_2d","test_out2")
test_result=compare_files("out2","test_out2")
if test_result:
	print("TEST OK\n")
else:
	print("TEST FAILED !!!")
	sys.exit(1)

print("Test model #3:")
cmd="ester 1d -M 10 -Xc 0.5 -i test_model1 -ndomains 16 -npts 20" 
cmd=cmd+" -noplot -tol 1e-8 -maxit 100 -o test_model3"
exec_cmd(cmd)
cmd="ester output test_model3"
exec_cmd(cmd,"template_1d","test_out3")
test_result=compare_files("out3","test_out3")
if test_result:
	print("TEST OK\n")
else:
	print("TEST FAILED !!!")
	sys.exit(1)

print("Test model #4:")
cmd="ester 2d -i test_model3 -Omega_bk 0.3" 
cmd=cmd+" -p 2d.par -noplot -tol 1e-8 -maxit 10 -o test_model4"
exec_cmd(cmd)
cmd="ester output test_model4"
exec_cmd(cmd,"template_2d","test_out4")
test_result=compare_files("out4","test_out4")
if test_result:
	print("TEST OK\n")
else:
	print("TEST FAILED !!!")
	sys.exit(1)

os.remove("test_model1")
os.remove("test_model2")
os.remove("test_model3")
os.remove("test_model4")
os.remove("test_out1")
os.remove("test_out2")
os.remove("test_out3")
os.remove("test_out4")

print("\n---------- All tests OK ---------------\n");
