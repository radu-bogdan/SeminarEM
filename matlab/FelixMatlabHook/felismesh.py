import sys
import numpy as npy
import time

if len(sys.argv) > 1:
    order = int(sys.argv[1])
    refine = int(sys.argv[2])
else:
    order = 1
    refine = 0


import numpy as np
from ngsolve import *
from netgen.occ import *
# from netgen.webgui import Draw as DrawGeo

print('Generating geometry:'); tm = time.monotonic()
def drawMagnet1(k):
    m1xnew = m1[0]*cos(k*pi/4) -m1[1]*sin(k*pi/4)
    m1ynew = m1[0]*sin(k*pi/4) +m1[1]*cos(k*pi/4)
    m1new = (m1xnew,m1ynew,0)

    m2xnew = m2[0]*cos(k*pi/4) -m2[1]*sin(k*pi/4)
    m2ynew = m2[0]*sin(k*pi/4) +m2[1]*cos(k*pi/4)
    m2new = (m2xnew,m2ynew,0)

    m3xnew = m3[0]*cos(k*pi/4) -m3[1]*sin(k*pi/4)
    m3ynew = m3[0]*sin(k*pi/4) +m3[1]*cos(k*pi/4)
    m3new = (m3xnew,m3ynew,0)

    m4xnew = m4[0]*cos(k*pi/4) -m4[1]*sin(k*pi/4)
    m4ynew = m4[0]*sin(k*pi/4) +m4[1]*cos(k*pi/4)
    m4new = (m4xnew,m4ynew,0)

    a5xnew = a5[0]*cos(k*pi/4) -a5[1]*sin(k*pi/4)
    a5ynew = a5[0]*sin(k*pi/4) +a5[1]*cos(k*pi/4)
    a5new = (a5xnew,a5ynew,0)

    a6xnew = a6[0]*cos(k*pi/4) -a6[1]*sin(k*pi/4)
    a6ynew = a6[0]*sin(k*pi/4) +a6[1]*cos(k*pi/4)
    a6new = (a6xnew,a6ynew,0)

    a7xnew = a7[0]*cos(k*pi/4) -a7[1]*sin(k*pi/4)
    a7ynew = a7[0]*sin(k*pi/4) +a7[1]*cos(k*pi/4)
    a7new = (a7xnew,a7ynew,0)

    a8xnew = a8[0]*cos(k*pi/4) -a8[1]*sin(k*pi/4)
    a8ynew = a8[0]*sin(k*pi/4) +a8[1]*cos(k*pi/4)
    a8new = (a8xnew,a8ynew,0)

    #Draw magnet
    seg1 = Segment(m1new,m2new);
    seg2 = Segment(m2new,m3new);
    seg3 = Segment(m3new,m4new);
    seg4 = Segment(m4new,m1new);
    magnet1 = Face(Wire([seg1,seg2,seg3,seg4]))
    #Draw air around magnet
    air_seg1 = Segment(m1new,a5new)
    air_seg2 = Segment(a5new,a6new)
    air_seg3 = Segment(a6new,m2new)
    air_seg4 = Segment(m2new,m1new)
    air_magnet1_1 = Face(Wire([air_seg1,air_seg2,air_seg3,air_seg4]))
    air_seg5 = Segment(m4new,m3new)
    air_seg6 = Segment(m3new,a7new)
    air_seg7 = Segment(a7new,a8new)
    air_seg8 = Segment(a8new,m4new)
    air_magnet1_2 = Face(Wire([air_seg5,air_seg6,air_seg7,air_seg8]))

    return (magnet1,air_magnet1_1,air_magnet1_2)

def drawMagnet2(k):
    m5xnew = m5[0]*cos(k*pi/4) -m5[1]*sin(k*pi/4)
    m5ynew = m5[0]*sin(k*pi/4) +m5[1]*cos(k*pi/4)
    m5new = (m5xnew,m5ynew,0)

    m6xnew = m6[0]*cos(k*pi/4) -m6[1]*sin(k*pi/4)
    m6ynew = m6[0]*sin(k*pi/4) +m6[1]*cos(k*pi/4)
    m6new = (m6xnew,m6ynew,0)

    m7xnew = m7[0]*cos(k*pi/4) -m7[1]*sin(k*pi/4)
    m7ynew = m7[0]*sin(k*pi/4) +m7[1]*cos(k*pi/4)
    m7new = (m7xnew,m7ynew,0)

    m8xnew = m8[0]*cos(k*pi/4) -m8[1]*sin(k*pi/4)
    m8ynew = m8[0]*sin(k*pi/4) +m8[1]*cos(k*pi/4)
    m8new = (m8xnew,m8ynew,0)

    a1xnew = a1[0]*cos(k*pi/4) -a1[1]*sin(k*pi/4)
    a1ynew = a1[0]*sin(k*pi/4) +a1[1]*cos(k*pi/4)
    a1new = (a1xnew,a1ynew,0)

    a2xnew = a2[0]*cos(k*pi/4) -a2[1]*sin(k*pi/4)
    a2ynew = a2[0]*sin(k*pi/4) +a2[1]*cos(k*pi/4)
    a2new = (a2xnew,a2ynew,0)

    a3xnew = a3[0]*cos(k*pi/4) -a3[1]*sin(k*pi/4)
    a3ynew = a3[0]*sin(k*pi/4) +a3[1]*cos(k*pi/4)
    a3new = (a3xnew,a3ynew,0)

    a4xnew = a4[0]*cos(k*pi/4) -a4[1]*sin(k*pi/4)
    a4ynew = a4[0]*sin(k*pi/4) +a4[1]*cos(k*pi/4)
    a4new = (a4xnew,a4ynew,0)

    #Draw magnet
    seg1 = Segment(m5new,m6new);
    seg2 = Segment(m6new,m7new);
    seg3 = Segment(m7new,m8new);
    seg4 = Segment(m8new,m5new);
    magnet2 = Face(Wire([seg1,seg2,seg3,seg4]))
    air_seg1 = Segment(m5new,a3new)
    air_seg2 = Segment(a3new,a4new)
    air_seg3 = Segment(a4new,m6new)
    air_seg4 = Segment(m6new,m5new)
    air_magnet2_1 = Face(Wire([air_seg1,air_seg2,air_seg3,air_seg4]))
    air_seg5 = Segment(m8new,m7new)
    air_seg6 = Segment(m7new,a2new)
    air_seg7 = Segment(a2new,a1new)
    air_seg8 = Segment(a1new,m8new)
    air_magnet2_2 = Face(Wire([air_seg5,air_seg6,air_seg7,air_seg8]))

    return (magnet2,air_magnet2_1,air_magnet2_2)

def drawStatorNut(k):
    s1xnew = s1[0]*cos(k*pi/24) -s1[1]*sin(k*pi/24)
    s1ynew = s1[0]*sin(k*pi/24) +s1[1]*cos(k*pi/24)
    s1new = (s1xnew,s1ynew,0)

    s2xnew = s2[0]*cos(k*pi/24) -s2[1]*sin(k*pi/24)
    s2ynew = s2[0]*sin(k*pi/24) +s2[1]*cos(k*pi/24)
    s2new = (s2xnew,s2ynew,0)

    s3xnew = s3[0]*cos(k*pi/24) -s3[1]*sin(k*pi/24)
    s3ynew = s3[0]*sin(k*pi/24) +s3[1]*cos(k*pi/24)
    s3new = (s3xnew,s3ynew,0)

    s4xnew = s4[0]*cos(k*pi/24) -s4[1]*sin(k*pi/24)
    s4ynew = s4[0]*sin(k*pi/24) +s4[1]*cos(k*pi/24)
    s4new = (s4xnew,s4ynew,0)

    s5xnew = s5[0]*cos(k*pi/24) -s5[1]*sin(k*pi/24)
    s5ynew = s5[0]*sin(k*pi/24) +s5[1]*cos(k*pi/24)
    s5new = (s5xnew,s5ynew,0)

    s6xnew = s6[0]*cos(k*pi/24) -s6[1]*sin(k*pi/24)
    s6ynew = s6[0]*sin(k*pi/24) +s6[1]*cos(k*pi/24)
    s6new = (s6xnew,s6ynew,0)

    s7xnew = s7[0]*cos(k*pi/24) -s7[1]*sin(k*pi/24)
    s7ynew = s7[0]*sin(k*pi/24) +s7[1]*cos(k*pi/24)
    s7new = (s7xnew,s7ynew,0)

    s8xnew = s8[0]*cos(k*pi/24) -s8[1]*sin(k*pi/24)
    s8ynew = s8[0]*sin(k*pi/24) +s8[1]*cos(k*pi/24)
    s8new = (s8xnew,s8ynew,0)

    #Draw stator coil
    seg1 = Segment(s2new,s3new)
    seg2 = Segment(s3new,s4new)
    seg3 = Segment(s4new,s5new)
    seg4 = Segment(s5new,s6new)
    seg5 = Segment(s6new,s7new)
    seg6 = Segment(s7new,s2new)
    stator_coil = Face(Wire([seg1,seg2,seg3,seg4,seg5,seg6]))
    #Draw air nut in the stator
    air_seg1 = Segment(s1new,s2new)
    air_seg2 = Segment(s2new,s7new)
    air_seg3 = Segment(s7new,s8new)
    air_seg4 = Segment(s8new,s1new)
    stator_air = Face(Wire([air_seg1,air_seg2,air_seg3,air_seg4]))

    stator_air = stator_air-(stator_air*air_gap_stator)
    return (stator_coil,stator_air)

def f48(s):
    return (s-1)%48

orign = (0,0);
#inner radius rotor
r1 = 26.5*10**(-3);
#outer radius rotor
r2 = 78.63225*10**(-3);
#sliding mesh rotor
r4 = 78.8354999*10**(-3);
#sliding mesh stator
r6 = 79.03874999*10**(-3);
#inner radius stator
r7 = 79.242*10**(-3);
#outer radius stator
r8 = 116*10**(-3)

#Points for magnet1 and air around magnet1
m1 = (69.23112999*10**(-3),7.535512*10**(-3),0)
m2 = (74.828958945*10**(-3),10.830092744*10**(-3),0)
m3 = (66.13621099700001*10**(-3),25.599935335*10**(-3),0)
m4 = (60.53713*10**(-3),22.30748*10**(-3),0)
a5 = (69.75636*10**(-3),5.749913*10**(-3),0)
a6 = (75.06735*10**(-3),3.810523*10**(-3),0)
a7 = (65.3506200*10**(-3),26.51379*10**(-3),0)
a8 = (59.942145092*10**(-3),24.083661604*10**(-3),0)
#Points for magnet2 and air around magnet2
m5 = (58.579985516*10**(-3), 27.032444757*10**(-3),0)
m6 = (64.867251151*10**(-3),28.663475405*10**(-3),0)
m7 = (60.570096319*10**(-3),45.254032279*10**(-3),0)
m8 = (54.282213127*10**(-3),43.625389857*10**(-3),0)
a1 = (53.39099766*10**(-3),45.259392713*10**(-3),0)
a2 = (55.775078884*10**(-3),50.386185578*10**(-3),0)
a3 = (59.41521771*10**(-3),25.355776837*10**(-3),0)
a4 = (65.12210917100001*10**(-3),27.707477175*10**(-3),0)
#Points for Stator Nut and air in the stator
s1 = (79.04329892000*10**(-3),3.9538335974*10**(-3),0)
s2 = (80.143057128*10**(-3),4.0037794254*10**(-3),0)
s3 = (80.387321219*10**(-3),2.965459706*10**(-3),0)
s4 = (98.78501315600001*10**(-3),3.9007973292*10**(-3),0)
s5 = (98.44904989600001*10**(-3),9.026606148400001*10**(-3),0)
s6 = (80.086666706*10**(-3),7.5525611543*10**(-3),0)
s7 = (79.980020247*10**(-3),6.4912415424*10**(-3),0)
s8 = (78.88229587*10**(-3),6.4102654448*10**(-3),0)

domains = []

# h_max = 0.001

# h_air_gap = 0.25*h_max
# h_air_magnets = h_max
# h_coils = h_max
# h_stator_air = h_max
# h_magnets = h_max
# h_stator_iron = h_max
# h_rotor_iron = h_max
# h_shaft_iron = h_max


h_max = 0.005

h_air_gap = 0.05*h_max
h_air_magnets = h_max
h_coils = h_max
h_stator_air = h_max
h_magnets = h_max
h_stator_iron = h_max
h_rotor_iron = h_max
h_shaft_iron = h_max

rotor_inner  = Circle(orign,r=r1).Face()
rotor_outer  = Circle(orign,r=r2).Face()
sliding_inner  = Circle(orign,r=r4).Face()
sliding_outer  = Circle(orign,r=r6).Face()
stator_inner = Circle(orign,r=r7).Face()
stator_outer = Circle(orign,r=r8).Face()

rotor_inner.edges[0].name = "rotor_inner"
rotor_outer.edges[0].name = "rotor_outer"
stator_inner.edges[0].name = "stator_inner"
stator_outer.edges[0].name = "stator_outer"

rotor_iron = rotor_outer - rotor_inner

air_gap_stator = stator_inner - sliding_outer
air_gap = sliding_outer - sliding_inner
air_gap_rotor = sliding_inner - rotor_outer

stator_iron = stator_outer - stator_inner

string_coils = ""
domain_name_coil = {0: "coil1", 1: "coil2", 2: "coil3", 3:"coil4", 4:"coil5",
                    5:"coil6", 6:"coil7", 7:"coil8", 8: "coil9", 9: "coil10", 10: "coil11", 11:"coil12",
                    12:"coil13",13: "coil14", 14: "coil15", 15: "coil16", 16:"coil17", 17:"coil18",
                    18: "coil19", 19: "coil20", 20: "coil21", 21:"coil22", 22:"coil23",
                    23: "coil24", 24: "coil25", 25: "coil26", 26:"coil27", 27:"coil28",
                    28: "coil29", 29: "coil30", 30: "coil31", 31:"coil32", 32:"coil33",
                    33: "coil34", 34: "coil35", 35: "coil36", 36:"coil37", 37:"coil38",
                    38: "coil39", 39: "coil40", 40: "coil41", 41:"coil42", 42:"coil43",
                    43: "coil44", 44: "coil45", 45: "coil46", 46:"coil47", 47:"coil48"}

for k in range(48):#48
    (stator_coil,stator_air) = drawStatorNut(k)

    stator_coil.faces.name = domain_name_coil[k]
    stator_air.faces.name = "air"

    stator_iron -= stator_coil
    stator_iron -= stator_air

    domains.append(stator_coil)
    domains.append(stator_air)
    string_coils += domain_name_coil[k] + "|"

domains_magnets = [];
domains_air_magnets = [];

domain_name_magnet1 = {0: "magnet1", 1: "magnet3", 2: "magnet5", 3:"magnet7", 4:"magnet9",
                    5:"magnet11", 6:"magnet13", 7:"magnet15"}
domain_name_magnet2 = {0: "magnet2", 1: "magnet4", 2: "magnet6", 3:"magnet8", 4:"magnet10",
                    5:"magnet12", 6:"magnet14", 7:"magnet16"}

for k in range(8):#8
    (magnet1,air_magnet1_1,air_magnet1_2) = drawMagnet1(k)
    (magnet2,air_magnet2_1,air_magnet2_2) = drawMagnet2(k)

    magnet1.faces.name = domain_name_magnet1[k]
    magnet1.faces.maxh = h_magnets
    magnet1.edges[0].name = "magnets_interface"
    magnet1.edges[1].name = "magnets_interface"
    magnet1.edges[2].name = "magnets_interface"
    magnet1.edges[3].name = "magnets_interface"

    magnet2.faces.name = domain_name_magnet2[k]
    magnet2.faces.maxh = h_magnets
    magnet2.edges[0].name = "magnets_interface"
    magnet2.edges[1].name = "magnets_interface"
    magnet2.edges[2].name = "magnets_interface"
    magnet2.edges[3].name = "magnets_interface"

    air_magnet1_1.faces.name = "rotor_air"
    air_magnet1_2.faces.name = "rotor_air"
    air_magnet2_1.faces.name = "rotor_air"
    air_magnet2_2.faces.name = "rotor_air"

    air_magnet1_1.faces.maxh = h_air_magnets
    air_magnet1_2.faces.maxh = h_air_magnets
    air_magnet2_1.faces.maxh = h_air_magnets
    air_magnet2_2.faces.maxh = h_air_magnets

    rotor_iron -= magnet1
    rotor_iron -= air_magnet1_1;
    rotor_iron -= air_magnet1_2;
    rotor_iron -= magnet2
    rotor_iron -= air_magnet2_1;
    rotor_iron -= air_magnet2_2;

    domains.append(magnet1)
    domains.append(magnet2)
    domains.append(air_magnet1_1)
    domains.append(air_magnet1_2)
    domains.append(air_magnet2_1)
    domains.append(air_magnet2_2)

stator_iron.faces.name = "stator_iron"
stator_iron.faces.maxh = h_stator_iron

air_gap_stator.faces.maxh = h_air_gap
air_gap_stator.faces.name = "air_gap_stator"

air_gap.faces.maxh = h_air_gap
air_gap.faces.name = "air_gap"

air_gap_rotor.faces.maxh = h_air_gap
air_gap_rotor.faces.name = "air_gap_rotor"

rotor_iron.faces.name = "rotor_iron"
rotor_iron.faces.maxh = h_rotor_iron

shaft_iron = rotor_inner
shaft_iron.faces.name = "shaft_iron"
shaft_iron.faces.maxh = h_shaft_iron

domains.append(shaft_iron)
domains.append(rotor_iron)
domains.append(air_gap_stator)
domains.append(air_gap)
domains.append(air_gap_rotor)
domains.append(stator_iron)

geo = Glue(domains)

geoOCC = OCCGeometry(geo, dim=2)

print('Done! took',time.monotonic()-tm)
tm = time.monotonic(); print('Generating mesh:')

geoOCCmesh = geoOCC.GenerateMesh()

if order==2:
    geoOCCmesh.SecondOrder()

if refine==1:
    geoOCCmesh.Refine()
    
    
print('Done! took',time.monotonic()-tm)
tm = time.monotonic(); print('Saving to mat...')

#####################################################################################################################
npoints = geoOCCmesh.Elements2D().NumPy()['np'].max()

p = geoOCCmesh.Coordinates()
t = npy.c_[geoOCCmesh.Elements2D().NumPy()['nodes'].astype(npy.uint64)[:,:npoints],
           geoOCCmesh.Elements2D().NumPy()['index'].astype(npy.uint64)]-1

e = npy.c_[geoOCCmesh.Elements1D().NumPy()['nodes'].astype(npy.uint64)[:,:((npoints+1)//2)],
           geoOCCmesh.Elements1D().NumPy()['index'].astype(npy.uint64)]-1

max_bc_index = geoOCCmesh.Elements1D().NumPy()['index'].astype(npy.uint64).max()
max_rg_index = geoOCCmesh.Elements2D().NumPy()['index'].astype(npy.uint64).max()

regions_1d_np = []
for i in range(max_bc_index):
    regions_1d_np += [geoOCCmesh.GetBCName(i)]

regions_2d_np = []
for i in range(max_rg_index):
    regions_2d_np += [geoOCCmesh.GetMaterial(i+1)]

identifications = npy.array(geoOCCmesh.GetIdentifications())
#####################################################################################################################



#####################################################################################################################
Mperp_mag1 = np.array([-0.507223091788922, 0.861814791678634])
Mperp_mag2 = np.array([-0.250741225095427, 0.968054150364350])
Mperp_mag3 = (-1)*np.array([-0.968055971101187, 0.250734195544481])
Mperp_mag4 = (-1)*np.array([-0.861818474866413, 0.507216833690415])
Mperp_mag5 = np.array([-0.861814791678634, -0.507223091788922])
Mperp_mag6 = np.array([-0.968054150364350, -0.250741225095427])
Mperp_mag7 = (-1)*np.array([-0.250734195544481, -0.968055971101187])
Mperp_mag8 = (-1)*np.array([-0.507216833690415, -0.861818474866413])
Mperp_mag9 = np.array([0.507223091788922, -0.861814791678634])
Mperp_mag10 = np.array([0.250741225095427, -0.968054150364350])
Mperp_mag11 = (-1)*np.array([0.968055971101187, -0.250734195544481])
Mperp_mag12 = (-1)*np.array([0.861818474866413, -0.507216833690415])
Mperp_mag13 = np.array([0.861814791678634, 0.507223091788922])
Mperp_mag14 = np.array([0.968054150364350, 0.250741225095427])
Mperp_mag15 = (-1)*np.array([0.250734195544481, 0.968055971101187])
Mperp_mag16 = (-1)*np.array([0.507216833690415, 0.861818474866413])

Mperp_mag = np.c_[Mperp_mag1,Mperp_mag2, Mperp_mag3, Mperp_mag4, Mperp_mag5, Mperp_mag6, Mperp_mag7, Mperp_mag8,
                  Mperp_mag9,Mperp_mag10,Mperp_mag11,Mperp_mag12,Mperp_mag13,Mperp_mag14,Mperp_mag15,Mperp_mag16]

m = np.c_[Mperp_mag[1,:],-Mperp_mag[0,:]].T
nu0 = 10**7/(4*np.pi)
m = m*nu0*1.158095238095238
#####################################################################################################################




#####################################################################################################################
offset = 0
polepairs  = 4
gamma_correction_model = -30.0
gamma = 40.0
gamma_correction_timestep = -1
phi0 = (gamma + gamma_correction_model + gamma_correction_timestep * polepairs) * np.pi/180.0


area_coils_UPlus  = np.r_[f48(offset+1) , f48(offset+2)]
area_coils_VMinus = np.r_[f48(offset+3) , f48(offset+4)]
area_coils_WPlus  = np.r_[f48(offset+5) , f48(offset+6)]
area_coils_UMinus = np.r_[f48(offset+7) , f48(offset+8)]
area_coils_VPlus  = np.r_[f48(offset+9) , f48(offset+10)]
area_coils_WMinus = np.r_[f48(offset+11), f48(offset+12)]

for k in range(1,4):
    area_coils_UPlus = np.r_[area_coils_UPlus, f48(k*12+offset+1)]
    area_coils_UPlus = np.r_[area_coils_UPlus, f48(k*12+offset+2)]
    
    area_coils_VMinus = np.r_[area_coils_VMinus, f48(k*12+offset+3)]
    area_coils_VMinus = np.r_[area_coils_VMinus, f48(k*12+offset+4)]
    
    area_coils_WPlus = np.r_[area_coils_WPlus, f48(k*12+offset+5)]
    area_coils_WPlus = np.r_[area_coils_WPlus, f48(k*12+offset+6)]
    
    area_coils_UMinus = np.r_[area_coils_UMinus, f48(k*12+offset+7)]
    area_coils_UMinus = np.r_[area_coils_UMinus, f48(k*12+offset+8)]
    
    area_coils_VPlus = np.r_[area_coils_VPlus, f48(k*12+offset+9)]
    area_coils_VPlus = np.r_[area_coils_VPlus, f48(k*12+offset+10)]
    
    area_coils_WMinus = np.r_[area_coils_WMinus, f48(k*12+offset+11)]
    area_coils_WMinus = np.r_[area_coils_WMinus, f48(k*12+offset+12)]
    
    
I0peak = 1555.63491861 ### *1.5
phase_shift_I1 = 0.0
phase_shift_I2 = 2/3*np.pi#4/3*pi
phase_shift_I3 = 4/3*np.pi#2/3*pi

I1c = I0peak * np.sin(phi0 + phase_shift_I1)
I2c = (-1)* I0peak * np.sin(phi0 + phase_shift_I2)
I3c = I0peak * np.sin(phi0 + phase_shift_I3)

areaOfOneCoil = 0.00018053718538758062

UPlus  =  I1c* 2.75 / areaOfOneCoil
VMinus =  I2c* 2.75 / areaOfOneCoil

WPlus  =  I3c* 2.75 / areaOfOneCoil
UMinus = -I1c* 2.75 / areaOfOneCoil

VPlus  = -I2c* 2.75 / areaOfOneCoil
WMinus = -I3c* 2.75 / areaOfOneCoil

j3 = np.zeros(48)
j3[area_coils_UPlus]  = UPlus
j3[area_coils_VMinus] = VMinus
j3[area_coils_WPlus]  = WPlus
j3[area_coils_UMinus] = UMinus
j3[area_coils_VPlus]  = VPlus
j3[area_coils_WMinus] = WMinus
#####################################################################################################################



#####################################################################################################################
import scipy.io
scipy.io.savemat('motor.mat', {"t" : t.T+1,
                               "e" : e.T+1,
                               "p" : p.T,
                               "regions_2d" : regions_2d_np,
                               "regions_1d" : regions_1d_np,
                               "m" : m,
                               "j3" : j3},
                 do_compression=True)
#####################################################################################################################

print('Done! took',time.monotonic()-tm)
