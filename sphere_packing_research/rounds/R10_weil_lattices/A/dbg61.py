import sys, json
from weil import *
p = 61
c = Ctx(); c.p = p; c.n = 120
c.B = np.load(os.path.join(HERE, 'B_p61_L0.npy')); c.E = np.load(os.path.join(HERE, 'E_p61_L0.npy'))
sp = math.sqrt(p)
c.Gp = (c.B + c.E / sp) / 2; c.Gm = (c.B - c.E / sp) / 2
c.D = 13725.0; c.theta = (p - 1) ** 2 / (4 * p)
F = np.load(os.path.join(HERE, 'F_p61_C3.npy'))
print('F shape', F.shape, 'max|F|', np.abs(F).max(), 'rank', np.linalg.matrix_rank(F.astype(float)))
o = grid_enum(c, F, label='p=61 C_3-fixed RECHECK')
x = o['best'][1]
xo = np.array([int(z) for z in x], dtype=object)
b = int(xo.dot(c.B.astype(object)).dot(xo)); e = int(xo.dot(c.E.astype(object)).dot(xo))
print('best x max|x|', max(abs(int(z)) for z in x), 'b', b, 'e', e, 'N', p * b * b - e * e, 'THR', (p - 1) ** 2 * 13725, 'Gamma exact', (p*b*b-e*e)/(4*p*13725))
# list all distinct gammas found near the threshold
json.dump(dict(x=[int(z) for z in x], b=b, e=e), open(os.path.join(HERE, 'dbg61_best.json'), 'w'))
