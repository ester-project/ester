#!/usr/bin/env python

import sys
import signal
import os
from subprocess import call

ester_root=ESTER_ROOT


# Let the child process handle SIGINT
def signal_handler(signal,frame):
	pass
signal.signal(signal.SIGINT, signal_handler);

command='help'
if len(sys.argv)>=2:
	command=sys.argv[1]

if command=='1d':
	cmd="star1d"
	call([ester_root+"/bin/"+cmd]+sys.argv[2:])

elif command=='2d':
	cmd="star2d"
	call([ester_root+"/bin/"+cmd]+sys.argv[2:])

elif command=='output':
	cmd="gen_output"
	call([ester_root+"/bin/"+cmd]+sys.argv[2:])
	
elif command=='info':
	cmd="ester_info"
	call([ester_root+"/bin/"+cmd]+sys.argv[2:])
	
elif command=='version':
	cmd="version"
	call([ester_root+"/bin/"+cmd]+sys.argv[2:])
	
elif command=='help':
	if len(sys.argv)<3:
		call(['cat',ester_root+"/doc/help/help"])
	else:
		subcmd=sys.argv[2]
		status=call(['cat',ester_root+"/doc/help/"+subcmd],stderr=open(os.devnull))
		if status:
			print("No help on '"+subcmd+"'")
else:
	print("Unknown command '"+command+"'")
