import numpy as np
import math
import calc

def calc_distance(self, point1, point2, rcutoff):
    """
        Calculates the distance between two points. If the the distance on one axis exceeds the cutoff distance, only that axes distance is returned.
    """

    deltaX = abs(point2[0] - point1[0])
    if deltaX > rcutoff:  return deltaX
    deltaY = abs(point2[1] - point1[1])
    if deltaY > rcutoff:  return deltaY
    deltaZ = abs(point2[2] - point1[2])
    if deltaZ > rcutoff:  return deltaZ

    return math.sqrt(np.power(deltaX,2) + deltaY**2 + deltaZ**2)

def calc_distance_f(self, point1, point2, rcutoff):
        calc.distance(point1, point2, rcutoff)
