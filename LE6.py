import numpy as np
import matplotlib.pyplot as plt

def f_system(x,y):
    if abs(x*y)<=1:
        f1=np.arcsin(x*y)
    else:
        f1=np.sign(x*y)*np.pi/2
    f2=np.exp(x+2*y-3)-1
    return f1, f2

def euler_system(f, x0, y0, b, n):
    h = (b - 0) / n
    x = [x0]
    y = [y0]
    for i in range(n):
        dx, dy = f(x[-1], y[-1])
        x_i = x[-1] + h * dx
        y_i = y[-1] + h * dy
        x.append(x_i)
        y.append(y_i)
    return x, y

def add_arrow(ax, x, y, dx, dy, color='blue', scale=0.1):
    norm = np.sqrt(dx**2 + dy**2)
    if norm > 0:
        ax.arrow(x, y, dx / norm * scale, dy / norm * scale,head_width=scale * 0.3, head_length=scale * 0.5, fc=color, ec=color, alpha=0.8)

def linearized_M(x, y):
    u = x - 0
    v = y - 1.5
    f1 = 1.5 * u
    f2 = u + 2 * v
    return f1, f2

def linearized_N(x, y):
    u = x - 3
    v = y - 0
    f1 = 3 * v
    f2 = u + 2 * v
    return f1, f2

plt.figure(figsize=(10, 8))

b_forward = 2.0
b_backward = -2.0
n = 150

x_vals_M1 = np.linspace(-1.0, 1.0, 100)
y_dir1 = -2 * x_vals_M1 + 1.5
plt.plot(x_vals_M1, y_dir1, 'r-', linewidth=2, alpha=0.8)

x_arrow = 0.6
y_arrow = -2 * x_arrow + 1.5
dx, dy = linearized_M(x_arrow, y_arrow)
add_arrow(plt.gca(), x_arrow, y_arrow, dx, dy, color='red', scale=0.15)

x_arrow2 = -0.6
y_arrow2 = -2 * x_arrow2 + 1.5
dx2, dy2 = linearized_M(x_arrow2, y_arrow2)
add_arrow(plt.gca(), x_arrow2, y_arrow2, dx2, dy2, color='red', scale=0.15)

plt.plot([0, 0], [-1, 3.5], 'b-', linewidth=2, alpha=0.8)
y_arrow_up = 2.8
dx_up, dy_up = linearized_M(0, y_arrow_up)
add_arrow(plt.gca(), 0, y_arrow_up, dx_up, dy_up, color='blue', scale=0.15)
y_arrow_down = 0.8
dx_down, dy_down = linearized_M(0, y_arrow_down)
add_arrow(plt.gca(), 0, y_arrow_down, dx_down, dy_down, color='blue', scale=0.15)

uzel_starts = [
    (0.1, 1.6), (0.1, 1.4),
    (-0.1, 1.6), (-0.1, 1.4),
    (0.2, 1.3), (-0.2, 1.7),
    (0.3, 1.2), (-0.3, 1.8),
]

for x0, y0 in uzel_starts:
    x_vals, y_vals = euler_system(linearized_M, x0, y0, b_forward, n)
    plt.plot(x_vals, y_vals, 'green', alpha=0.7, linewidth=1.5)
    x_vals_back, y_vals_back = euler_system(linearized_M, x0, y0, b_backward, n)
    plt.plot(x_vals_back, y_vals_back, 'green', alpha=0.7, linewidth=1.5)

    if len(x_vals) > 30:
        for idx in [15, 35]:
            if idx < len(x_vals):
                dx, dy = linearized_M(x_vals[idx], y_vals[idx])
    add_arrow(plt.gca(), x_vals[idx], y_vals[idx], dx, dy, color='darkgreen', scale=0.1)

x_vals1 = np.linspace(0.5, 5.0, 100)
y_vals1 = -x_vals1 / 3 + 1
plt.plot(x_vals1, y_vals1, 'r-', linewidth=2, alpha=0.8)

x_arrow_left = 2.2
y_arrow_left = -x_arrow_left / 3 + 1
dx_left, dy_left = linearized_N(x_arrow_left, y_arrow_left)
add_arrow(plt.gca(), x_arrow_left, y_arrow_left, dx_left, dy_left, color='red', scale=0.15)


x_arrow_right = 3.8
y_arrow_right = -x_arrow_right / 3 + 1
dx_right, dy_right = linearized_N(x_arrow_right, y_arrow_right)
add_arrow(plt.gca(), x_arrow_right, y_arrow_right, dx_right, dy_right, color='red', scale=0.15)

x_vals2 = np.linspace(1.5, 4.5, 100)
y_vals2 = x_vals2 - 3
plt.plot(x_vals2, y_vals2, 'b-', linewidth=2, alpha=0.8)

x_arrow_right_up = 3.8
y_arrow_right_up = x_arrow_right_up - 3
dx_right_up, dy_right_up = linearized_N(x_arrow_right_up, y_arrow_right_up)
add_arrow(plt.gca(), x_arrow_right_up, y_arrow_right_up, dx_right_up, dy_right_up, color='blue', scale=0.15)

x_arrow_left_down = 2.4
y_arrow_left_down = x_arrow_left_down - 3
dx_left_down, dy_left_down = linearized_N(x_arrow_left_down, y_arrow_left_down)
add_arrow(plt.gca(), x_arrow_left_down, y_arrow_left_down, dx_left_down, dy_left_down, color='blue', scale=0.15)

sedlo_starts = [
    (2.9, 0.3), (3.1, 0.3),
    (2.9, -0.3), (3.1, -0.3),

]

for x0, y0 in sedlo_starts:
    x_vals, y_vals = euler_system(linearized_N, x0, y0, b_forward, n//2)
    plt.plot(x_vals, y_vals, 'purple', alpha=0.7, linewidth=1.5)
    x_vals_back, y_vals_back = euler_system(linearized_N, x0, y0, b_backward, n//2)
    plt.plot(x_vals_back, y_vals_back, 'purple', alpha=0.7, linewidth=1.5)

    if len(x_vals) > 30:
        for idx in [20, 40]:
            if idx < len(x_vals):
                dx, dy = linearized_N(x_vals[idx], y_vals[idx])
    idx=15
    add_arrow(plt.gca(), x_vals[idx], y_vals[idx], dx, dy, color='darkviolet', scale=0.08)

plt.plot(0, 1.5, 'ro', markersize=12, markeredgewidth=2,markerfacecolor='red', markeredgecolor='darkred', label='M(0, 1.5) - неустойчивый узел')
plt.plot(3, 0, 'go', markersize=12, markeredgewidth=2,markerfacecolor='lime', markeredgecolor='darkgreen', label='N(3, 0) - седло')
plt.axhline(0, color='black', linewidth=0.5, alpha=0.3)
plt.axvline(0, color='black', linewidth=0.5, alpha=0.3)
plt.xlabel("x", fontsize=12)
plt.ylabel("y", fontsize=12)
plt.title("Фазовые портреты линеаризованных систем", fontsize=14)
plt.grid(True, alpha=0.2)
plt.legend(fontsize=10, loc='upper right')
plt.axis('equal')
plt.xlim(-1, 3.5)
plt.ylim(-1, 3.5)
plt.tight_layout()
plt.show()

