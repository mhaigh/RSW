import multiprocessing as mp
import numpy as np
import sys
import os
import time

#========================================================================

pe = 2;

def foo(pi):
	for i in range(0,2):
		return(pi,i);	

def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

if len(sys.argv) > 1:
	ncpus = int(sys.argv[1])
else:
	ncpus = 1;

def main():
	p = mp.Pool(pe);
	for result in p.imap(foo,[1,2]):
		x=result
		

	


if __name__ == '__main__':
	main();

print(x)
