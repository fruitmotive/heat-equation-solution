import numpy as np
import matplotlib.pyplot as plt

xs, ys, ts = 0, 0, 0
xe, ye, te = 10, 5, 0.2
nx, ny, nt = 100, 100, 100
x, y, t = np.linspace(xs, xe, nx), np.linspace(ys, ye, ny), np.linspace(ts, te, nt)
hx, hy, tau = x[1] - x[0], y[1] - y[0], t[1] - t[0]
gx, gy = tau / hx**2, tau / hy**2
u = np.zeros((nx, ny, 2*nt + 1))
x, y = np.meshgrid(x, y)
u[:, :, 0] = np.sin(np.pi*x) * np.cos(2*np.pi*y)

def foox(i1, i2, j):
    return 0.5*gy*(u[i1, i2-1, j-1] + u[i1, i2+1, j-1]) + (1 - gy)*u[i1, i2, j-1]

def fooy(i1, i2, j):
    return 0.5*gx*(u[i1-1, i2, j-1] + u[i1+1, i2, j-1]) + (1 - gx)*u[i1, i2, j-1]

def fx(i2, j):
    d, s = np.zeros(nx), np.zeros(nx)
    d[1] = 0
    s[1] = 0
    A = 0.5*gx
    B = 1 + gx
    C = A
    for m in range(1, nx - 1):
        foom = -foox(m, i2, j)
        d[m+1] = C/(B - A*d[m])
        s[m+1] = (foom - A*s[m])/(A*d[m] - B)
    u[nx-1, i2, j] = 0
    for m in range(nx - 1, 0, -1):
        u[m-1, i2, j] = d[m] * u[m, i2, j] + s[m]

def fy(i1, j):
    d, s = np.zeros(ny), np.zeros(ny)
    d[1] = 1
    s[1] = 0
    A = 0.5*gy
    B = 1 + gy
    C = A
    for m in range(1, ny - 1):
        foom = -fooy(i1, m, j)
        d[m+1] = C/(B - A*d[m])
        s[m+1] = (foom - A*s[m])/(A*d[m] - B)
    u[i1, ny-1, j] = s[-1] / (1 - d[-1])
    for m in range(ny - 1, 0, -1):
        u[i1, m-1, j] = d[m] * u[i1, m, j] + s[m]

for j in range(1, 2*nt, 2):
    for i2 in range(1, ny-1):
        fx(i2, j)
    for i1 in range(1, nx-1):
        fy(i1, j+1)
    
ax = plt.axes(projection="3d")
ax.plot_surface(x, y, u[:, :, 150], cmap='plasma')
ax.set_zlim(-2, 2)
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("graphresh150.png", dpi=400)

