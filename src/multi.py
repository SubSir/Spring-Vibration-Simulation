import ctypes
import os

# 获取当前进程 ID
pid = os.getpid()
# 设置进程优先级为高优先级
handle = ctypes.windll.kernel32.OpenProcess(0x00000010, False, pid)
ctypes.windll.kernel32.SetPriorityClass(handle, 0x00000080)  # HIGH_PRIORITY_CLASS

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sympy as sp

# 定义常量
n = 12  # 振子数量
length = 1.0
length0 = 0.5
mass = 1.0
kk = 10.0
m_last = 2.0  # 最后一个振子的质量

k = kk * n  # 弹簧劲度系数
m = mass / (n - 1)  # 每个振子的质量（除了最后一个）
l0 = length0 / n  # 弹簧的原长

# 定义符号变量
t = sp.symbols("t")
q = [sp.Function(f"q{i}")(t) for i in range(1, n + 1)]  # 位移 q1(t), q2(t), ..., qn(t)
q_dot = [sp.diff(qi, t) for qi in q]  # 速度 q1'(t), q2'(t), ..., qn'(t)

# 定义势能和动能
V = 0.5 * k * (q[0] - l0) ** 2  # 第一个振子的势能
T = 0.0

for i in range(n):
    T += 0.5 * (m if i < n - 1 else m_last) * q_dot[i] ** 2  # 动能
    if i > 0:
        V += 0.5 * k * (q[i] - q[i - 1] - l0) ** 2  # 中间弹簧的势能

# 定义拉格朗日量
L = T - V

# 计算拉格朗日方程
lagrange_eqs = [sp.diff(sp.diff(L, q_dot[i]), t) - sp.diff(L, q[i]) for i in range(n)]

# 解出二阶导数 q1_ddot, q2_ddot, ..., qn_ddot
q_ddot = [sp.solve(eq, qi.diff(t, t))[0] for eq, qi in zip(lagrange_eqs, q)]


# 定义状态向量 [q1, q1_dot, q2, q2_dot, ..., qn, qn_dot]
def ode_func(t, y):
    y_dict = {q[i]: y[2 * i] for i in range(n)}
    for i in range(n):
        y_dict[q_dot[i]] = y[2 * i + 1]
    y_ddot = [q_ddot[i].subs(y_dict) for i in range(n)]
    ans = []
    for i in range(n):
        ans.append(y[2 * i + 1])
        ans.append(y_ddot[i])
    return ans


# 定义初始条件
initial_conditions = []
for i in range(n):
    initial_conditions.append(length / n * (i + 1))  # 初始位置
    initial_conditions.append(0.0)  # 初始速度

# 时间跨度
t_span = (0, 10)
t_eval = np.linspace(t_span[0], t_span[1], 100000)

# 解ODE
solution = solve_ivp(ode_func, t_span, initial_conditions, t_eval=t_eval)

# 绘制结果
plt.figure(figsize=(10, 5))
for i in range(n):
    plt.plot(solution.t, solution.y[2 * i], label=f"Position (q{i+1})")

time = 0
times = 0
for i in range(1, len(solution.t) - 1):
    if (
        solution.y[2 * (n - 1)][i] > solution.y[2 * (n - 1)][i - 1]
        and solution.y[2 * (n - 1)][i] > solution.y[2 * (n - 1)][i + 1]
    ):
        times += 1
        time = solution.t[i]

print("TIME: ", time / times)
plt.xlabel("Time (s)")
plt.ylabel("Position")
plt.title("Spring-Mass System Simulation with n Groups")
plt.legend()
plt.grid(True)
plt.show()
