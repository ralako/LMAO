import numpy as np
import matplotlib.pyplot as plt
import math
import ROS
from celluloid import Camera


ROS.colormesh("50x50",10,T=200,fps=50,smooth=False,uloz=True)    # rozmery, cas,...
#ROS.colormesh_anim("40x40","Acka",200,range(500),50,uloz=True)

#ROS.colormesh_BD("50x50",10,T=200,fps=50,uloz=True)        # rozmery, cas,...
#ROS.colormesh_BD_anim("50x50",T=200,fps=50,uloz=True)

#ROS.colormesh_T("50x50",10,r=1,fps=50,smooth=False,uloz=True)  # zavislost na T
#ROS.colormesh_T_anim("50x50",range(30),fps=1)
