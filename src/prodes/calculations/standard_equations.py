import math
import numpy as np


def distance(a, b):
    """calculates the distance between point a and b in 3d space"""

    dist = math.sqrt((a.x - b.x)**2 + (a.y - b.y)**2 + (a.z - b.z)**2)
    return dist


def direction(a, b):
    """calculates the vector from point a to b in 3d space"""
    x = (b.x - a.x)
    y = (b.y - a.y)
    z = (b.z - a.z)
    return x, y, z


def find_middel(a, b):
    """calculates the coordinates at the middel between point a
       and b in 3d space"""
    x = (a.x + b.x)/2
    y = (a.y + b.y)/2
    z = (a.z + b.z)/2
    return x, y, z


def trimean(values):
    """Calculates the trimean of a set of values"""

    q2 = np.median(values)
    q1 = np.median([x for x in values if x < q2])
    q3 = np.median([x for x in values if x > q2])
    return (q2*2 + q1 + q3)/4

def pos_charge(pka,ph):
    """calculate the charge of a group which can be protonated"""

    return 1/(1+10**(ph-pka))


def neg_charge(pka,ph):
    """calculate the charge of a group which can be deprotonated"""

    return -1/(1+10**(pka-ph))

def angle_d2(vector1, vector2):
    """Calculates the angle between two 3d vectors in radians"""
    vec1_norm, vec2_norm = np.linalg.norm(vector1), np.linalg.norm(vector2)
    cos_theta = np.dot(vector1, vector2)/(vec1_norm*vec2_norm)
    return np.arccos(cos_theta)


def rotation_x(theta):
  return np.matrix([[ 1, 0           , 0           ],
                   [ 0, math.cos(theta),-math.sin(theta)],
                   [ 0, math.sin(theta), math.cos(theta)]])
  
def rotation_y(theta):
  return np.matrix([[ math.cos(theta), 0, math.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-math.sin(theta), 0, math.cos(theta)]])
  
def rotation_z(theta):
  return np.matrix([[ math.cos(theta), -math.sin(theta), 0 ],
                   [ math.sin(theta), math.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])