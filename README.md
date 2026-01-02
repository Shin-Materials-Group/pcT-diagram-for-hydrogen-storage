# pcT-diagram-for-hydrogen-storage (code)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import Normalize
def draw_pct_diagram_for_input_temp(temperature_list=None):
    Delta_H = 26956.183
    Delta_S = 102.754
    f_s = 0.575

    A_alpha = 5.229E-4
    gamma_alpha = 1.103
    V_alpha = -1.05E-2
    H_alpha = -11573.73

    A_beta = 3.572E-2
    gamma_beta = 0.546
    V_beta = -1.22E-5
    H_beta = -13619.08

    R = 8.314

    C_range_min, C_range_max, C_step = 0, 2.2, 0.005
    P_min, P_max = 0.01, 10000

    C_range = np.arange(C_range_min, C_range_max + C_step, C_step)
    P_guess_range = np.linspace(P_min, P_max, 100000)
    
    fig, ax = plt.subplots(figsize=(10,7))

    if temperature_list is None:
        print("온도를 입력하세요(°C). 입력 없이 바로 엔터 누르면 종료합니다.")
        temps_to_plot = []
        while True:
            temp_input = input("온도 (°C): ")
            if temp_input.strip() == '':
                break
            try:
                temp_val = float(temp_input)
                temps_to_plot.append(temp_val)
                
            except:
                print("숫자를 입력해 주세요.")
    else:
        temps_to_plot = temperature_list

    cmap = colormaps['viridis']
    norm = Normalize(vmin=0, vmax=len(temps_to_plot)-1)
    
    for idx, T_C in enumerate(temps_to_plot):
        T = T_C + 273.15
        P_ref = 1.0

        def P_plateau(C):
            return P_ref * np.exp((-Delta_H/(R*T)) + (Delta_S/R) + f_s * (C))

        def C_single_phase(P, A, gamma, V, H):
            return A * (P ** (gamma/2)) * np.exp(-gamma * V * P / (R * T)) * np.exp(-gamma * H / (R * T))

        def find_crossing_C(A, gamma, V, H):
            min_diff = float('inf')
            crossing_C = None
            for C in C_range:
                Peq = P_plateau(C)
                diff = np.abs(C_single_phase(P_guess_range, A, gamma, V, H) - C)
                P_phase = P_guess_range[np.argmin(diff)]
                diff_P = np.abs(P_phase - Peq)
                diff_C = np.abs(C_single_phase(P_phase, A, gamma, V, H) - C)
                total_diff = diff_P + diff_C
                if total_diff < min_diff:
                    min_diff = total_diff
                    crossing_C = C
            return crossing_C

        C_alpha_eq = find_crossing_C(A_alpha, gamma_alpha, V_alpha, H_alpha)
        C_beta_eq = find_crossing_C(A_beta, gamma_beta, V_beta, H_beta)

        P_final = []
        for C in C_range:
            P_eq = P_plateau(C)
            diff_alpha = np.abs(C_single_phase(P_guess_range, A_alpha, gamma_alpha, V_alpha, H_alpha) - C)
            P_alpha = P_guess_range[np.argmin(diff_alpha)]
            diff_beta = np.abs(C_single_phase(P_guess_range, A_beta, gamma_beta, V_beta, H_beta) - C)
            P_beta = P_guess_range[np.argmin(diff_beta)]
            if C < C_alpha_eq:
                P_final.append(min(P_alpha, P_eq))
            elif C_alpha_eq <= C <= C_beta_eq:
                P_final.append(P_eq)
            else:
                P_final.append(max(P_beta, P_eq))

        color = cmap(norm(idx))
        ax.semilogy(C_range, P_final, label=f'{T_C}°C', color=color)

    ax.set_xlim([C_range_min, C_range_max])
    ax.set_ylim([P_min, P_max])
    ax.set_xlabel('Concentration C (wt%)')
    ax.set_ylabel('Pressure P (bar)')
    ax.set_title('PCT Diagram (Accumulated Graph by Temperature)')
    ax.legend()
    ax.grid(True, which="both")
    plt.show()

Main contributor: Jiheon Yoo (Dankook Univ.)
