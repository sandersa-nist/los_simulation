# -----------------------------------------------------------------------------
# Name:   simulation.py     
# Purpose:    A simple line of site power simulator
# Authors:     aric.sanders@nist.gov
# Created:     03/19/205
# License:     NIST License
# -----------------------------------------------------------------------------
"""docstring
"""
# -----------------------------------------------------------------------------
# Standard Imports
import sys
import os

# -----------------------------------------------------------------------------
# Third Party Imports
sys.path.append(os.path.join(os.path.dirname( __file__ ),'..')) # if the repo is library/source/modules, one deep change to .
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy import ndimage


# -----------------------------------------------------------------------------
# Module Constants
# -----------------------------------------------------------------------------
# Module Functions
def simple_directional_gain(theta_array,theta1=-(np.pi/180)*30,theta2=(np.pi/180)*30,gain1=9,gain2=-100): 
    return np.where((theta_array>=theta1)&(theta_array<=theta2),gain1,gain2)

def omni(theta_array,gain=9):
    return gain*np.ones(theta_array.size)

def cos_three_halves(theta_array,gain=1):
    return gain+10*np.log10(np.cos(theta_array)**(3/2))

def ntia_very_high_gain_model_point(theta,gain =50 ):
    """From TM 09-461 """
    assert gain>48, "gain is too low to use this model, ntia_high_gain or ntia_medium_gain"
    theta = abs(theta)
    theta_m = 50*(.25*gain+7)**(.5)/10**(gain/20)
    theta_r = 27.466*10**(-.3*gain/10)
    theta_b = 48
    gain_out  = 0
    if 0<=theta<np.pi/180 *theta_m:
        gain_out = gain - 4*10**(-4)*(10**(gain/10))*(theta*180/np.pi)**2
    elif np.pi/180 *theta_m<=theta<np.pi/180 *theta_r:
        gain_out = .75*gain - 7
    elif np.pi/180 *theta_r<=theta<np.pi/180 *theta_b:
        gain_out = 29-25*np.log10(theta*180/np.pi)
    else:
        gain_out=-13
    return gain_out
ntia_very_high_gain_model = np.vectorize(ntia_very_high_gain_model_point,excluded=set(["gain"]))

def ntia_high_gain_model_point(theta,gain =30):
    """From TM 09-461 """
    assert 22<=gain<48, "gain is too low or too high to use this model, try ntia_very_high_gain or ntia_medium_gain"
    theta = abs(theta)
    theta_m = 50*(.25*gain+7)**(.5)/10**(gain/20)
    theta_r = 250/10**(gain/20)
    theta_b = 48
    gain_out  = 0
    if 0<=theta<np.pi/180 *theta_m:
        gain_out = gain - 4*10**(-4)*(10**(gain/10))*(theta*180/np.pi)**2
    elif np.pi/180 *theta_m<=theta<np.pi/180 *theta_r:
        gain_out = .75*gain - 7
    elif np.pi/180 *theta_r<=theta<np.pi/180 *theta_b:
        gain_out = 29-25*np.log10(theta*180/np.pi)
    else:
        gain_out=-13
    return gain_out
ntia_high_gain_model = np.vectorize(ntia_high_gain_model_point,excluded=set(["gain"]))

def ntia_medium_gain_model_point(theta,gain =20 ):
    """From TM 09-461 """
    assert 10<=gain<22, "gain is too low or too high to use this model, try ntia_very_high_gain or ntia_medium_gain"
    theta = abs(theta)
    theta_m = 50*(.25*gain+7)**(.5)/10**(gain/20)
    theta_r = 250/10**(gain/20)
    theta_b = 131.8257 * 10**(-gain/50)
    gain_out  = 0
    if 0<=theta<np.pi/180 *theta_m:
        gain_out = gain - 4*10**(-4)*(10**(gain/10))*(theta*180/np.pi)**2
    elif np.pi/180 *theta_m<=theta<np.pi/180 *theta_r:
        gain_out = .75*gain - 7
    elif np.pi/180 *theta_r<=theta<np.pi/180 *theta_b:
        gain_out = 53-gain/2-25*np.log10(theta*180/np.pi)
    else:
        gain_out = 0
    return gain_out
ntia_medium_gain_model = np.vectorize(ntia_medium_gain_model_point,excluded=set(["gain"]))

def calculate_horizon(height):
    """Calculates the horizon distance in m given the height in meters"""
    return 3.56972*10**3*np.sqrt(height)


def calculate_relative_angle(x1,y1,x2,y2):
        mag_a = np.sqrt(x2**2+y2**2)
        mag_b = np.sqrt(x1**2+y1**2)
        #print(f"mag_a:{mag_a}, mag_b:{mag_b}")
        #print((x1*x2+y1*y2)/(mag_a*mag_b))
        inverse_cos =  np.arccos((x1*x2+y1*y2)/(mag_a*mag_b))
        v1 = [x1,y1,0]
        v2= [x2,y2,0]
        cross_product =np.cross(v1,v2)
        if cross_product[2]<0:
            angle = -inverse_cos
        else:
            angle = inverse_cos
        return angle 
        
vcal = np.vectorize(calculate_relative_angle,excluded=set(["x1","y1"]))
def node_distance(node1,node2):
    return np.sqrt((node1.location[0]-node2.location[0])**2+(node1.location[1]-node2.location[1])**2)

def node_to_node_power(node1,node2,wavelength = 299792458/3.75e9):
    """Returns the loss in dB between 2 nodes"""
    distance = node_distance(node1,node2)
    rx_angle = node1.calculate_relative_angle(node2.location[0]-node1.location[0],node2.location[1]-node1.location[1])
    tx_angle = node2.calculate_relative_angle(node1.location[0]-node2.location[1],node1.location[1]-node2.location[1])
    gain_rx = node1.antenna_pattern(rx_angle)
    gain_tx = node2.antenna_pattern(tx_angle)
    power_tx = node2.power
    power_rx = power_tx + gain_rx+ gain_tx + 20*np.log10(wavelength/(4*np.pi*distance))
    if isinstance(power_rx,np.ndarray):
        power_rx = power_rx[0]
    return power_rx

def node_to_node_loss(node1,node2,wavelength = 299792458/3.75e9):
    """Returns the loss in dB between 2 nodes"""
    distance = node_distance(node1,node2)
    rx_angle = node1.calculate_relative_angle(node2.location[0]-node1.location[0],node2.location[1]-node1.location[1])
    tx_angle = node2.calculate_relative_angle(node1.location[0]-node2.location[1],node1.location[1]-node2.location[1])
    gain_rx = node1.antenna_pattern(rx_angle)
    gain_tx = node2.antenna_pattern(tx_angle)
    loss = gain_rx+ gain_tx + 20*np.log10(wavelength/(4*np.pi*distance))
    return loss

# -----------------------------------------------------------------------------
# Module Classes
class Node():
    def __init__(self,direction,location,antenna_pattern=simple_directional_gain,id=None):
        self.direction = direction
        self.location = location
        self.antenna_pattern = antenna_pattern
        self.id = id

    def calculate_relative_angle(self,x,y):
        return vcal(self.direction[0],self.direction[1],x,y)
    
    def plot_antenna_pattern(self,**options):
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        theta = np.linspace(-np.pi, np.pi, 100)
        relative_theta = self.calculate_relative_angle(np.cos(theta),np.sin(theta))
        r = self.antenna_pattern(relative_theta)
        theta_direction = np.arctan2(self.direction[1],self.direction[0])
        ax.plot(theta_direction+relative_theta,r,color="b",marker="D",linestyle="solid")
        ax.set_title(f"Antenna Pattern Gain", va='bottom')
        ax.grid(True)



class TxNode(Node):
    def __init__(self,direction,location,power = 1,antenna_pattern=simple_directional_gain,signal=1,id=None):
        super().__init__(direction=direction,location=location,antenna_pattern=antenna_pattern,id=id)
        self.signal = signal
        self.power = power
        
    def calculate_relative_angle(self,x,y):
        return vcal(self.direction[0],self.direction[1],x,y)

        
    
class RxNode(Node):
    def __init__(self,direction,location,antenna_pattern=simple_directional_gain,id=None):
        super().__init__(direction=direction,location=location,antenna_pattern=antenna_pattern,id=id)
# -----------------------------------------------------------------------------
# Module Scripts

def plot_antenna_functions(antenna_functions=[omni,simple_directional_gain,
                                              cos_three_halves,
                                              ntia_very_high_gain_model,ntia_high_gain_model,
                                              ntia_medium_gain_model]):
    """plots all antenna functions in antenna_functions"""
    figure,ax = plt.subplots(subplot_kw={'projection': 'polar'})
    theta = np.linspace(-1*np.pi,1*np.pi,1000)
    for func in antenna_functions:
        try:
            ax.plot(theta,func(theta),label = func.__name__)
        except:
            pass
    ax.legend()

# -----------------------------------------------------------------------------
# Module Runner
if __name__=="__main__":
    plot_antenna_functions()