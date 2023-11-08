import numpy as np
from rdkit import Chem
from math import sqrt, atan2, pi, ceil
from xml.etree import ElementTree


def get_torsion_library(xml):

    tree = ElementTree.parse(xml)
    xml_tree = tree.getroot()

    patterns = {}

    rank = 0

    for c in xml_tree.findall("hierarchyClass"):
        cname = c.get("name")

        patterns[cname] = []

        for r in c.iter("torsionRule"):
            rmethod = r.get("method")
            rsmarts = r.get("smarts")
            rsmarts_mol = Chem.MolFromSmarts(rsmarts)

            if rmethod == "exact":

                renergy = []
                rci_lower = []
                rci_upper = []

                for b in r.find("histogram_converted").findall("bin"):
                    benergy = float(b.get("energy"))
                    bci_lower = float(b.get("lower"))
                    bci_upper = float(b.get("upper"))

                    renergy.append(benergy)
                    rci_lower.append(bci_lower)
                    rci_upper.append(bci_upper)

                patterns[cname].append({"method":rmethod,
                                        "smarts":rsmarts, 
                                        "smarts_mol": rsmarts_mol,
                                        "energy":renergy, 
                                        "ci_lower":rci_lower, 
                                        "ci_upper":rci_upper})

            else:

                rtheta_0 = []
                rtolerance2 = []
                rbeta_1 = []
                rbeta_2 = []

                for a in r.find("angleList").findall("angle"):
                    atheta_0 = float(a.get("theta_0"))
                    atolerance2 = float(a.get("tolerance2"))
                    abeta_1 = float(a.get("beta_1"))
                    abeta_2 = float(a.get("beta_2"))

                    rtheta_0.append(atheta_0)
                    rtolerance2.append(atolerance2)
                    rbeta_1.append(abeta_1)
                    rbeta_2.append(abeta_2)

                patterns[cname].append({"method":rmethod,
                                        "smarts":rsmarts, 
                                        "smarts_mol": rsmarts_mol,
                                        "theta_0": rtheta_0,
                                        "tolerance2":rtolerance2,
                                        "beta_1":rbeta_1,
                                        "beta_2":rbeta_2})

            if cname == "GG":
                patterns[cname][-1]["rank"] = rank + 9999
            else:
                patterns[cname][-1]["rank"] = rank    

            rank += 1

    return patterns


def ang_diff(theta_1, theta_2):
    # (-180,180] -> [0, 360)
    if theta_1 < 0:
        theta_1 += 360
    if theta_2 < 0:
        theta_2 += 360
    del_theta = (theta_1 - theta_2) % 360 #Angular difference
    # [0, 360) -> (-180, 180]
    if del_theta > 180:
        del_theta -= 360
    return(del_theta)
    # Test this works: ang_diff(0, 179), ang_diff(0, -179),
    # and ang_diff(-179, 179)


def calculate_score(theta, scores):

    bin_idx = ceil(theta / 10) + 17
    score = (scores[bin_idx]-scores[(bin_idx+35)%36])/10.0*(theta-(bin_idx-17)*10)+scores[bin_idx]
    return score

def unit(a):
    # The argument should be a NumPy array with 1 axis
    return(a / sqrt(np.dot(a,a))) #Scales a by its norm


def dihedral(a_1, a_2, a_3, a_4):
    # The arguments should all be NumPy arrays with 1 axis and length 3
    # These atoms should be in order, with a_2 and a_3 defining the bond
    # of interest

    # The 3 displacement vectors:
    b_1 = a_2 - a_1
    b_2 = a_3 - a_2
    b_3 = a_4 - a_3

    n_1 = unit(np.cross(b_1, b_2))
    n_2 = unit(np.cross(b_2, b_3))

    # Imagine the first atom (a_1) is above the middle bond (from a_2 to a_3),
    # so that b_1 points downward. Then n_1 points out of the page

    m = unit(np.cross(n_1, b_2))
    # I moved the normalization to be after the cross product. Moving
    # the normalization should not change the end result because the cross
    # product commutes with scalar multiplication and ||n_1|| = 1

    # Looking down b_2, we can consider n_1 to be the x-axis and
    # m to be the y-axis. Then the dihedral angle is the angle that
    # n_2 makes with the x-axis when projected into this plane. Since dihedral
    # angles are measured going clockwise, we need to negate the angle
    # that we get back from atan
    x = np.dot(n_1, n_2) #Project n_2 onto n_1
    y = np.dot(m, n_2) #Project n_2 onto m
    return(-atan2(y,x) * 180 / pi) #Return the angle in degrees


def process_match(theta, pattern):

    if pattern["method"] == "exact": 
            
        energy = calculate_score(theta, pattern["energy"])
        lower_energy = calculate_score(theta, pattern["ci_lower"])
        upper_energy = calculate_score(theta, pattern["ci_upper"])

    else:
       
        energy = None

        for i, theta_0 in enumerate(pattern["theta_0"]):
            delta = ang_diff(theta, theta_0)

            if abs(delta) <= pattern["tolerance2"][i]:
                beta_1 = pattern["beta_1"][i]
                beta_2 = pattern["beta_2"][i]

                energy = beta_1 * (delta ** 2) + beta_2 * (delta ** 4)
                break

        lower_energy = energy
        upper_energy = energy

    return {"energy":energy, "ci_lower":lower_energy, "ci_upper":upper_energy}