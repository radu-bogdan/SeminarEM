import numpy as np

motor_npz = np.load('motor.npz', allow_pickle=True)

p = motor_npz['p']
e = motor_npz['e']
t = motor_npz['t']
regions_2d = motor_npz['regions_2d']
regions_1d = motor_npz['regions_1d']
m = motor_npz['m']
j3 = motor_npz['j3']

