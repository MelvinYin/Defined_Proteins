###############################################################################
# geometry

# A set of geometry functions for manipulating pdb files.
###############################################################################

from math import sqrt, cos, sin, acos, pi
import numpy as np
import scipy.spatial
import math

def crossProduct(u, v):
    """
    Calculates the cross product of two 3d vectors (as 1-d arrays).
    """
    return np.cross(u, v)  # possible need to convert from np.array

def dotProduct(u, v):
    """
    Calculates the dot product between two vectors.
    """
    return np.dot(u, v)

def findAngle(u, v):
    """
    Calculates the angle (degrees) between two vectors (as 1-d arrays) using
    dot product.
    """
    mag_u = np.linalg.norm(u)
    mag_v = np.linalg.norm(v)

    return np.rad2deg(math.acos(np.dot(u, v) / (mag_u * mag_v)))

def calcDihedrals(prevCO, currN, currCA, currCO, nextN):
    """
    Calculates phi and psi angles for an individual residue.
    Requires coord tuple of each atom.
    """
    prevCO = np.array(prevCO)
    currN = np.array(currN)
    currCA = np.array(currCA)
    currCO = np.array(currCO)
    nextN = np.array(nextN)

    # Set CA coordinates to origin
    A = prevCO - currCA
    B = currN - currCA
    C = currCO - currCA
    D = nextN - currCA

    # Calculate necessary cross products (define vectors normal to planes)
    V1 = crossProduct(A, B)
    V2 = crossProduct(C, B)
    V3 = crossProduct(C, D)

    # Determine scalar angle between normal vectors
    phi = findAngle(V1, V2)  # angle between pCO-CA X N-CA, and CO-CA X N-CA
    if dotProduct(A, V2) > 0:  # positive means pCO-CA and CO-CA X N-CA in same
        # direction
        phi = -phi

    psi = findAngle(V2, V3)  # angle between CO-CA X N-CA and CO-CA X nN-CA
    # positive means nN-CA and CO-CA X N-CA in same direction
    if dotProduct(D, V2) < 0:
        psi = -psi
    return phi, psi

def genRotMatrix(axis, theta):
    """
    Generate a rotation matrix for rotation of theta about axis.
    """

    matrix = np.zeros((3, 3), dtype=float)

    axis_length = sqrt((axis[0] ** 2 + axis[1] ** 2 + axis[2] ** 2))
    xNorm = axis[0] / axis_length
    yNorm = axis[1] / axis_length
    zNorm = axis[2] / axis_length

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    one_costheta = 1.0 - cos_theta

    matrix[0][0] = cos_theta + xNorm * xNorm * one_costheta
    matrix[0][1] = xNorm * yNorm * one_costheta - zNorm * sin_theta
    matrix[0][2] = xNorm * zNorm * one_costheta + yNorm * sin_theta
    matrix[1][0] = xNorm * yNorm * one_costheta + zNorm * sin_theta
    matrix[1][1] = cos_theta + yNorm * yNorm * one_costheta
    matrix[1][2] = yNorm * zNorm * one_costheta - xNorm * sin_theta
    matrix[2][0] = xNorm * zNorm * one_costheta - yNorm * sin_theta
    matrix[2][1] = yNorm * zNorm * one_costheta + xNorm * sin_theta
    matrix[2][2] = cos_theta + zNorm * zNorm * one_costheta

    return matrix

###############################################################################
# Unused, deprecated

# def get_length(vector):
#     return np.linalg.norm(vector)
#
# def dist(c1,c2):
#     """
#     Calculate the distance between two coordinates in 3d space.
#     """
#     return scipy.spatial.distance.euclidean(c1, c2)
#
# def dist_sq(c1,c2):
#     """
#     Calculate the squared distance between two coordinates in 3d space.
#     """
#     return scipy.spatial.distance.cdist([c1], [c2])[0] ** 2

# def calcDistances(coord):
#     """
#     Calculate all distances in coord.
#     """
#     return scipy.spatial.distance.cdist(coord, coord)
#     # possible need to convert from np.array

# def genRotMatrix_original(axis, theta):
#     """
#     Return the rotation matrix associated with clockwise rotation about
#     the given axis by theta degree.
#     """
#     axis = np.asarray(axis)
#     axis = axis/math.sqrt(np.dot(axis, axis))
#     a = math.cos(theta/2.0)
#     b, c, d = axis*math.sin(theta/2.0)
#     aa, bb, cc, dd = a*a, b*b, c*c, d*d
#     bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
#     return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
#                      [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
#                      [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

# def arbRot(matrix, axis, theta):
#     return np.dot(rotation_matrix(axis,theta), matrix)
#
# def arbRotCoord(coord,axis,theta):
#     """
#     Rotate all vectors in coord about an arbitray axis by theta.
#     """
#     return [[np.dot(rotation_matrix(axis,theta), c)] for c in coord]


# def calcGlyCbeta(Ncoord,CAcoord,COcoord):
#     """
#     Generates a beta carbon for a glycine using the coordinates for the amide
#     N, alpha C, and carboxyl C.
#     """
#
#     CA_CO_vector = []; CA_N_vector = []
#     for i in range(3):
#         CA_CO_vector.append(COcoord[i] - CAcoord[i])
#         CA_N_vector.append(Ncoord[i] - CAcoord[i])
#
#     rotation_amount = 240*(pi/180.)
#
#     rotated = arbRot(CA_CO_vector, CA_N_vector, rotation_amount)
#
#     CBeta = []
#     for i in range(3):
#         CBeta.append(rotated[i] + CAcoord[i])
#
#     return CBeta


# def calcHXT(C_coord,O_coord,OXT_coord):
#     """
#     Calculates the location of HXT using the location of C, O, and OXT.
#     (C-terminal hydrogen).
#     """
#     C_coord = np.array(C_coord)
#     O_coord = np.array(O_coord)
#     OXT_coord = np.array(OXT_coord)
#
#     O_C_dist = O_coord - C_coord
#     OXT_C_dist = OXT_coord - C_coord
#
#     C_O_OXT_vect = O_C_dist + OXT_C_dist
#     vect_normalised = C_O_OXT_vect/get_length(vector_3)
#
#     HXT_coord = OXT_coord + vect_normalised
#
#     return HXT_coord

# def calcHG(CB_coord,SG_coord):
#     """
#     Calculates the location of HG using the location of CB and SG.
#     (Hydrogen on free cysteines).
#     """
#     CB_coord = np.array(CB_coord)
#     SG_coord = np.array(SG_coord)
#
#     SG_CB_vect = SG_coord - CB_coord
#
#     vect_normalised = SG_CB_vect / get_length(SG_CB_vect)
#
#     HG_coord = SG_coord + 1.08 * vect_normalised
#
#     return HG_coord

# def calcHN(prevCO,prevO,currN):
#     """
#     Calculate the position of the amide hydrogen.
#     """
#     prevCO = np.array(prevCO)
#     prevO = np.array(prevO)
#     currN = np.array(currN)
#
#     CO_bond = prevO - prevCO
#     CN_bond = currN - prevCO
#
#     return prevCO + CO_bond + CN_bond

