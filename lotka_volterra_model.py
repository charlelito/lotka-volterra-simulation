# -*- coding: utf-8 -*-
"""
Model: Lotka-Volterra (Predator-Prey)
Generated from: Voltera.ipynb
Author: Charles Pimenta
charlelito@gmail.com
Exported at: 2025-10-31T12:40:21

Note:
- IPython magics ('%...') and shell lines ('!...') were commented out for .py validity.
- Each cell is separated by a banner indicating its original index and nearest preceding markdown heading.
"""


#==============================================================================
#  Cell 3 | Section: Implementação em Python
#==============================================================================
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


@dataclass
class LVParams:
    """Parâmetros do modelo de Lotka–Volterra."""
    alpha: float = 1.0   # taxa de crescimento das presas (sem predadores) [%]
    beta:  float = 0.1   # taxa de predação (acoplamento X–Y) [%]
    gamma: float = 1.5   # taxa de mortalidade natural dos predadores (sem presas) [%]
    delta: float = 0.075 # eficiência de conversão de presas em predadores ou taxa de reprodução dos predadores [%]

# Função diferencial
def lv_ode(t: float, z: np.ndarray, p: LVParams) -> list:
    """
    Campo vetorial do sistema de Lotka–Volterra.
    z = [X, Y], sendo X presas e Y predadores.
    """
    x, y = z
    dxdt = p.alpha * x - p.beta * x * y
    dydt = p.delta * x * y - p.gamma * y
    return [dxdt, dydt]

def simulate_lv(params: LVParams,
                z0=(10.0, 5.0),
                t_span=(0.0, 50.0),
                n_points=1000,
                rtol=1e-8, atol=1e-10):
    """
    Integra numericamente o sistema e retorna tempos e trajetórias (X(t), Y(t)).
    """
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(lambda t, z: lv_ode(t, z, params),
                    t_span=t_span, y0=np.array(z0, dtype=float),
                    t_eval=t_eval, rtol=rtol, atol=atol, vectorized=False)
    if not sol.success:
        raise RuntimeError(f"Falha na integração: {sol.message}")
    return sol.t, sol.y[0], sol.y[1]

def equilibrium_points(params: LVParams):
    """
    Retorna os dois pontos de equilíbrio:
    E0 = (0, 0)
    E1 = (gamma/delta, alpha/beta)
    """
    e0 = (0.0, 0.0)
    e1 = (params.gamma / params.delta, params.alpha / params.beta)
    return e0, e1

def plot_time_series(t, X, Y, params: LVParams, z0,filename='lv_time_series.png'):
    plt.figure(figsize=(9, 4.5))
    plt.plot(t, X, label='Presas (X)')
    plt.plot(t, Y, label='Predadores (Y)')
    plt.xlabel('Tempo')
    plt.ylabel('População')
    plt.title('Dinâmica Predador–Presa (Lotka–Volterra)')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=160)
    plt.show()

def plot_phase_plane(X, Y, params: LVParams, z0,filename='lv_phase_plane.png'):
    e0, e1 = equilibrium_points(params)

    plt.figure(figsize=(5.5, 5.5))
    plt.plot(X, Y, linewidth=2.0, label='Trajetória')
    plt.scatter([z0[0]], [z0[1]], marker='p', s=120,  color='red', label='Condição inicial em $t_0$')
    plt.scatter([e0[0]], [e0[1]], marker='x', s=80, label='Equilíbrio (0,0)')
    plt.scatter([e1[0]], [e1[1]], marker='o', s=60, label=r'Equilíbrio $(\gamma/\delta,\ \alpha/\beta)$')
    plt.xlabel('Presas (X)')
    plt.ylabel('Predadores (Y)')
    plt.title('Plano de Fase')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(filename, dpi=160)
    plt.show()

#==============================================================================
#  Cell 5 | Section: Gráficos: séries temporais e plano de fase
#==============================================================================
params = LVParams(alpha=1.0, beta=0.1, gamma=1.5, delta=0.075)
z0 = (10.0, 5.0)
t, X, Y = simulate_lv(params, z0=z0, t_span=(0, 25), n_points=1200)
plot_time_series(t, X, Y, params, z0,filename='lv_time_series_1.png')
plot_phase_plane(X, Y, params, z0,filename='lv_phase_plane_1.png')

#==============================================================================
#  Cell 7 | Section: Análise da sensibilidade a variação de parâmetros (variação de $\alpha$ e $\gamma$)
#==============================================================================
def sweep_params(base: LVParams, alphas = [0.8, 1.0, 1.2], gammas = [1.2, 1.5, 1.8], z0=(10.0, 5.0), t_span=(0, 40), n_points=800):
    # Variação de alpha (presas)
    plt.figure(figsize=(9, 4.5))
    for a in alphas:
        p = LVParams(alpha=a, beta=base.beta, gamma=base.gamma, delta=base.delta)
        t, X, Y = simulate_lv(p, z0=z0, t_span=t_span, n_points=n_points)
        plt.plot(t, X, label=fr'$\alpha={a:.1f}$')
    plt.xlabel('Tempo'); plt.ylabel('Presas (X)')
    plt.title('Sensibilidade: variação de α (presas)')
    plt.grid(True, alpha=0.3); plt.legend()
    plt.tight_layout(); plt.savefig('lv_sens_alpha.png', dpi=160); plt.show()
    
    # Plano de fase
    plt.figure(figsize=(5.5, 5.5))
    for a in alphas:
        p = LVParams(alpha=a, beta=base.beta, gamma=base.gamma, delta=base.delta)
        t, X, Y = simulate_lv(p, z0=z0, t_span=t_span, n_points=n_points)
        plt.plot(X, Y, linewidth=2.0, label=fr'Trajetória: $\alpha={a:.1f}$')
    plt.xlabel('Presas (X)')
    plt.ylabel('Predadores (Y)')
    plt.title('Plano de Fase')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('lv_phase_plane_var_alpha.png', dpi=160)
    plt.show()

    # Variação de gamma (predadores)
    plt.figure(figsize=(9, 4.5))
    for g in gammas:
        p = LVParams(alpha=base.alpha, beta=base.beta, gamma=g, delta=base.delta)
        t, X, Y = simulate_lv(p, z0=z0, t_span=t_span, n_points=n_points)
        plt.plot(t, Y, label=fr'$\gamma={g:.1f}$')
    plt.xlabel('Tempo'); plt.ylabel('Predadores (Y)')
    plt.title('Sensibilidade: variação de γ (predadores)')
    plt.grid(True, alpha=0.3); plt.legend()
    plt.tight_layout(); plt.savefig('lv_sens_gamma.png', dpi=160); plt.show()
    
    # Plano de fase
    plt.figure(figsize=(5.5, 5.5))
    for g in gammas:
        p = LVParams(alpha=base.alpha, beta=base.beta, gamma=g, delta=base.delta)
        t, X, Y = simulate_lv(p, z0=z0, t_span=t_span, n_points=n_points)
        plt.plot(X, Y, linewidth=2.0, label=fr'Trajetória: $\gamma={g:.1f}$')
    plt.xlabel('Presas (X)')
    plt.ylabel('Predadores (Y)')
    plt.title('Plano de Fase')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('lv_phase_plane_var_gamma.png', dpi=160)
    plt.show()
        

base = LVParams(alpha=1.0, beta=0.1, gamma=1.5, delta=0.075)
alphas = [0.8, 1.0, 1.2]
gammas = [1.2, 1.5, 1.8]
sweep_params(base,alphas,gammas)

#==============================================================================
#  Cell 9 | Section: Validação rápida dos equilíbrios
#==============================================================================
params = LVParams()
e0, e1 = equilibrium_points(params)
params = LVParams(alpha=1.0, beta=0.1, gamma=1.5, delta=0.075)

# Derivadas em cada equilíbrio (devem ser ~0):
for name, pt in [('E0', e0), ('E1', e1)]:
    dx, dy = lv_ode(0.0, np.array(pt, dtype=float), params)
    print(f"{name} = {pt} -> (dX/dt, dY/dt) = ({dx:.6f}, {dy:.6f})")
    
    t, X, Y = simulate_lv(params, z0=pt, t_span=(0, 25), n_points=1200)
    
    filename = f'lv_time_series_eq_{pt[0]:.0f}_{pt[1]:.0f}.png'
    plot_time_series(t, X, Y, params, pt,filename=filename)
    filename = f'lv_phase_plane_eq_{pt[0]:.0f}_{pt[1]:.0f}.png'
    plot_phase_plane(X, Y, params, pt,filename=filename)   

