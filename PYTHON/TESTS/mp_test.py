import multiprocessing as mp
import numpy as np
import sys
import os
import time

#========================================================================

def foo():
	for i in range(0,10):
		print(i);	

def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

if len(sys.argv) > 1:
	ncpus = int(sys.argv[1])
else:
	ncpus = 1;


if __name__ == '__main__':
	info('main line');
	#mp.set_start_method('spawn');
	p1 = mp.Process(target=foo)
	p1.start();
	p2 = mp.Process(target=foo)
	p2.start();
	p3 = mp.Process(target=foo)
	p3.start();
	p4 = mp.Process(target=foo)
	p4.start();
	p5 = mp.Process(target=foo)
	p5.start();
	p6 = mp.Process(target=foo)
	p6.start();
	p7 = mp.Process(target=foo)
	p7.start();
	p8 = mp.Process(target=foo)
	p8.start();


