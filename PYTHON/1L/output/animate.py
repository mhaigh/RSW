# animate.py
#=======================================================
#=======================================================

# File of input parameters for the 1L RSW plunger code

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from inputFile import *
#=======================================================

plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg';

filename = '/home/mike/Documents/GulfStream/RSW/DATA/1L/STOCH/u_nd.npy';
filename = '/home/mike/Documents/GulfStream/RSW/PYTHON/1L/u_nd.npy';

u_nd = np.load(filename);
wn = np.shape(u_nd)[2];
ulim = np.max(abs(u_nd[:,:,:]));
print(wn);
fig = plt.figure()
ax = fig.add_subplot(111);
#ax.set_xlim([x_grid.min(), x_grid.max()]);
#ax.set_ylim([y_grid.min(), y_grid.max()]);

# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
ims = []
for i in range(0,wn):
#for i in range(int(wn/2),wn):
#for i in range(wn-6,wn+6):
	ii = i%wn;
	print ii
	s = str(int(i/2)) + ' days';
	t = ax.annotate(s,(-0.4,-0.4));
	im = plt.imshow(u_nd[:,:,ii],animated=True,vmin=-ulim,vmax=ulim);
	#im = plt.pcolor(x_grid, y_grid, u_nd[:,:,ii], cmap='bwr', vmin=-ulim, vmax=ulim);
	##im.;
	ims.append([im,t])
	
ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True,repeat_delay=0)

ani.save('dynamic_images.mp4',metadata={'artist':'Guido'},writer='ffmpeg',fps=20)

plt.show()
