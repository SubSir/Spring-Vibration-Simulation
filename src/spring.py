import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sympy as sp

# 定义符号变量
t = sp.symbols("t")
q1 = sp.Function("q1")(t)  # 位移 q1(t)
q1_dot = sp.diff(q1, t)  # 速度 q1'(t)
q2 = sp.Function("q2")(t)  # 位移 q2(t)
q2_dot = sp.diff(q2, t)  # 速度 q2'(t)

# 定义常量
k = 10.0  # 弹簧劲度系数
m = 1.0  # 小球质量

# 定义势能和动能
V = 0.5 * k * (q1**2 + (q2 - q1) ** 2)  # 弹性势能
T = 0.5 * m * (q1_dot**2 + q2_dot**2)  # 动能

# 定义拉格朗日量
L = T - V

# 计算拉格朗日方程
lagrange_eq = [
    sp.diff(sp.diff(L, q1_dot), t) - sp.diff(L, q1),
    sp.diff(sp.diff(L, q2_dot), t) - sp.diff(L, q2),
]

# 解出二阶导数 q1_ddot 和 q2_ddot
q_ddot = [
    sp.solve(lagrange_eq[0], q1.diff(t, t))[0],
    sp.solve(lagrange_eq[1], q2.diff(t, t))[0],
]


# 定义状态向量 [q1, q1_dot, q2, q2_dot]
def ode_func(t, y):
    y_dict = {q1: y[0], q1_dot: y[1], q2: y[2], q2_dot: y[3]}
    q1_ddot_val = q_ddot[0].subs(y_dict)
    q2_ddot_val = q_ddot[1].subs(y_dict)
    return [y[1], q1_ddot_val, y[3], q2_ddot_val]


# 定义初始条件
initial_conditions = [1.0, 0.0, 0.0, 0.0]  # 初始位置和速度

# 时间跨度
t_span = (0, 10)
t_eval = np.linspace(t_span[0], t_span[1], 400)

# 解ODE
solution = solve_ivp(ode_func, t_span, initial_conditions, t_eval=t_eval)

# 绘制结果
plt.figure(figsize=(10, 5))
plt.plot(solution.t, solution.y[0], label="Position (q1)")
plt.plot(solution.t, solution.y[2], label="Position (q2)")
plt.xlabel("Time (s)")
plt.ylabel("Position")
plt.title("Spring-Mass System Simulation with Two Groups")
plt.legend()
plt.grid(True)
plt.show()
