from dataclasses import dataclass
from enum import Enum
import math as m, sys, multiprocessing, copy
from scipy import optimize


S_Fstat_min = 3.5
S_Fdyn_interval = (1.5, 2)
S_Hstat_min = 1.3
S_Hdyn_interval = (1.2, 1.5)

def involute(alpha):
    return m.tan(m.radians(alpha)) - m.radians(alpha)
def inverse_involute(alpha, anfangswert = 20):
    try:
        return float(optimize.newton(lambda x: involute(x) - alpha, anfangswert))
    except RuntimeError:
        assert(False)
   
Ritzel = 0
Rad = 1
indices = Ritzel, Rad

@dataclass
class Profil:
    alpha_n : int
    h_aP_s : float
    h_fP_s : float
    rho_fP_s : float
Normalprofil1 =     Profil(20, 1, 1.25, 0.375)
Normalprofil2 =     Profil(20, 1, 1.25, 0.375)
Protuberanzprofil = Profil(20, 1, 1.40, 0.400)

@dataclass
class Werkstoff:
    class Art(Enum):
        EH_Stahl = 1
        
    art : Art
    sigma_Hlim : int
    sigma_FE : int


def execute(P, n, b_d_1_verhaeltnis, verzahnungs_qualitaten, m_n, z_1, z_2, x_1, x_2, alpha_n, beta, k, sigma_Hlim, sigma_FE, K_A, K_S):
  assert(beta == 0)
  assert(k == 0)

  print("Konstanten")
  print(f"P = {P}")
  print(f"n = {n}")
  print(f"m_n = {m_n}")
  print(f"z_1 = {z_1}")
  print(f"z_2 = {z_2}")
  print(f"x_1 = {x_1}")
  print(f"x_2 = {x_2}")
  print(f"alpha_n = {alpha_n}")
  print(f"beta = {beta}")
  print(f"k = {k}")
  print(f"sigma_Hlim = {sigma_Hlim}")
  print(f"sigma_FE = {sigma_FE}")
  print(f"K_A = {K_A}")
  print(f"K_S = {K_S}")
  print()

  print("Geometrie und Kinetik")

  h_aP = m_n
  h_fP = 1.25*m_n
  rho_fP = 0.25*m_n
  print(f"h_aP = {h_aP}")
  print(f"h_fP = {h_fP}")
  print(f"rho_fP = {rho_fP}")

  def d(z):
    return m_n * z

  u = z_2 / z_1
  print(f"u = {u}")

  d_1 = d(z_1)
  d_2 = d(z_2)
  print(f"d_1 = {d_1}")
  print(f"d_2 = {d_2}")

  b = d_1 * b_d_1_verhaeltnis
  print(f"b = {b}")

  assert(all(q in (6, 7) for q in verzahnungs_qualitaten))
  assert(b / m_n <= 30) # Konstruktionsvorgaben, Tabelle 4

  def d_b(d):
    return d * m.cos(m.radians(alpha_n))
  d_b1 = d_b(d_1)
  d_b2 = d_b(d_2)
  print(f"d_b1 = {d_b1}")
  print(f"d_b2 = {d_b2}")

  def d_a(d, x):
    return d + 2 * (x*m_n+h_aP+k*m_n)
  d_a1 = d_a(d_1, x_1)
  d_a2 = d_a(d_2, x_2)
  print(f"d_a1 = {d_a1}")
  print(f"d_a2 = {d_a2}")

  def d_f(d, x):
    return d - 2 * (h_fP - x * m_n)
  d_f1 = d_f(d_1, x_1)
  d_f2 = d_f(d_2, x_2)
  print(f"d_f1 = {d_f1}")
  print(f"d_f2 = {d_f2}")

  h = h_aP + k*m_n+h_fP
  print(f"h = {h}")

  alpha_t = alpha_n
  alpha_wt = inverse_involute(involute(alpha_n) + 2 * (x_1 + x_2) / (z_1+z_2)*m.tan(m.radians(alpha_n)) )
  print(f"alpha_wt = {alpha_wt}")

  def d_w(d_b):
    return d_b / m.cos(m.radians(alpha_wt))
  d_w1 = d_w(d_b1)
  d_w2 = d_w(d_b2)
  print(f"d_w1 = {d_w1}")
  print(f"d_w2 = {d_w2}")

  a_w = 0.5 * (d_w1 + d_w2)
  print(f"a_w = {a_w}")

  # Profilueberdeckung

  p_n = m_n*m.pi
  print(f"p_n = {p_n}")

  p_t = m_n*m.pi/m.cos(m.radians(beta))
  print(f"p_t = {p_t}")

  epsilon_alpha = (m.sqrt(d_a1**2 - d_b1**2) + m.sqrt(d_a2**2 - d_b2**2) - (d_b1+d_b2)*m.tan(m.radians(alpha_wt))) / (2 * p_t * m.cos(m.radians(alpha_t)))
  print(f"epsilon_alpha = {epsilon_alpha}")

  epsilon_beta = 0
  print(f"epsilon_beta = {epsilon_beta}")

  epsilon_gamma = epsilon_alpha + epsilon_beta
  print(f"epsilon_gamma = {epsilon_gamma}")
  assert(epsilon_gamma > 1)

  # Unterschnitt

  h_FaP0 = h_fP / m_n
  print(f"h_FaP0 = {h_FaP0}")

  def z_min(x):
    return 2 * (h_FaP0 - x) / m.sin(m.radians(alpha_t))**2
  z_1min = z_min(x_1)
  z_2min = z_min(x_2)
  print(f"z_1min = {z_1min}")
  print(f"z_2min = {z_2min}")
  assert(z_1min <= z_1)
  assert(z_2min <= z_2)

  # Spitzwerden

  def gamma(z, x):
    return inverse_involute(m.pi / 2 / z + 2 * x * m.tan(m.radians(alpha_n)) / z + involute(alpha_t))
  gamma_1 = gamma(z_1, x_1)
  gamma_2 = gamma(z_2, x_2)
  print(f"gamma_1 = {gamma_1}")
  print(f"gamma_2 = {gamma_2}")

  def d_amax(z, gamma):
    return m_n * z * m.cos(m.radians(alpha_t)) / m.cos(m.radians(gamma))
  d_a1max = d_amax(z_1, gamma_1)
  d_a2max = d_amax(z_2, gamma_2)
  print(f"d_a1max = {d_a1max}")
  print(f"d_a2max = {d_a2max}")
  assert(d_a1max >= d_a1)
  assert(d_a2max >= d_a2)
  
  v = n / 60 * 2 * m.pi * d_1 / 2000
  print(f"v = {v}")

  F_t = 1000 * P / v
  print(f"F_t = {F_t}")
  
  T_1 = F_t * d_1 / 2000
  T_2 = F_t * d_2 / 2000
  print(f"T_1 = {T_1}")
  print(f"T_2 = {T_2}")


  # s_pr: Fussfreischnitt
  def Festigkeitsnachweis(verzahnungsqualitaet, s, s_pr, ist_ritzel):
    if verzahnungsqualitaet == 6:
      K_1 = 9.6
    elif verzahnungsqualitaet == 7:
      K_1 = 15.3
    else:
      raise NotImplementedError
    print(f"K_1 = {K_1}")
    if beta == 0:
      K_2 = 0.0193
    else:
      K_2 = 0.0087
    print(f"K_2 = {K_2}")

    temp1 = z_1 * v / 100 * m.sqrt(u**2 / (1 + u**2))
    assert(temp1 < 10)  # Glg 3.04

    temp2 = max(K_A * F_t / b, 100)
    K_V = 1 + ( K_1 / temp2 + K_2) * temp1
    print(f"K_V = {K_V}")

    assert(F_t / b * K_A >= 100)  # Abschnitt 3.4.1

    F_m = F_t * K_A * K_V
    print(f"F_m / b = {F_m / b}")

    assert(s == 0)
    A = 0.023

    f_sh = F_m / b * A * (b / d_1)**2
    print(f"f_sh = {f_sh}")

    f_ma = 0

    F_betax = abs(1.33 * f_sh)
    print(f"F_betax = {F_betax}")

    gamma_beta = 0.15 * F_betax
    print(f"gamma_beta = {gamma_beta}")

    F_betay = F_betax - gamma_beta
    print(f"F_betay = {F_betay}")

    c_gamma = 20

    K_Hbeta = 1 + c_gamma * F_betay / 2 / (F_m / b)
    if K_Hbeta > 2:
      K_Hbeta = m.sqrt(2 * c_gamma * F_betay / (F_m / b))
    print(f"K_Hbeta = {K_Hbeta}")

    K_Fbeta = K_Hbeta ** (1 / (1 + h / b + (h / b)**2))
    print(f"K_Fbeta = {K_Fbeta}")

    assert(verzahnungsqualitaet in (6, 7))
    K_Halpha = 1
    K_Falpha = 1
    print(f"K_Halpha = {K_Halpha}")
    print(f"K_Falpha = {K_Falpha}")

    # Gruebchentragfaehigkeit

    Z_H = m.sqrt((2 * m.cos(m.radians(beta)) * m.cos(m.radians(alpha_wt))) / (m.cos(m.radians(alpha_t))**2 * m.sin(m.radians(alpha_wt))))
    print(f"Z_H = {Z_H}")

    Z_E = 189.8
    print(f"Z_E = {Z_E}")

    E = 206000
    print(f"E = {E}")

    Z_epsilon = m.sqrt((4 - epsilon_alpha) / 3)
    print(f"Z_epsilon = {Z_epsilon}")

    Z_beta = m.sqrt(m.cos(m.radians(beta)))
    print(f"Z_beta = {Z_beta}")

    sigma_H0 = Z_H * Z_E * Z_epsilon * Z_beta * m.sqrt(F_t * (u + 1) / d_1 / b / u)
    print(f"sigma_H0 = {sigma_H0}")

    M1 = m.tan(m.radians(alpha_wt)) / m.sqrt((m.sqrt(d_a1**2/d_b1**2-1)-2*m.pi/z_1)*(m.sqrt(d_a2**2/d_b2**2-1)-(epsilon_alpha - 1)*2*m.pi/z_2))
    print(f"M1 = {M1}")

    M2 = m.tan(m.radians(alpha_wt)) / m.sqrt((m.sqrt(d_a2**2/d_b2**2-1)-2*m.pi/z_2)*(m.sqrt(d_a1**2/d_b1**2-1)-(epsilon_alpha - 1)*2*m.pi/z_1))
    print(f"M2 = {M2}")

    Z_B = 1
    if M1 > 1:
      Z_B = M1
    print(f"Z_B = {Z_B}")

    Z_D = 1
    if M2 > 1:
      Z_D = M2
    print(f"Z_D = {Z_D}")

    if ist_ritzel:  # ritzel
      sigma_Hdyn = Z_B * sigma_H0 * m.sqrt(K_A * K_V * K_Hbeta * K_Halpha)
      sigma_Hstat = Z_B * sigma_H0 * m.sqrt(K_S * K_V * K_Hbeta * K_Halpha)
    else: # rad
      sigma_Hdyn = Z_D * sigma_H0 * m.sqrt(K_A * K_V * K_Hbeta * K_Halpha)
      sigma_Hstat = Z_D * sigma_H0 * m.sqrt(K_S * K_V * K_Hbeta * K_Halpha)
    print(f"sigma_Hdyn = {sigma_Hdyn}")
    print(f"sigma_Hstat = {sigma_Hstat}")

    # zulaessige beanspruchung

    Z_NTstat = 1.6
    Z_NTdyn = 1
    print(f"Z_NTstat = {Z_NTstat}")
    print(f"Z_NTdyn = {Z_NTdyn}")

    if m_n <= 10:
      Z_X_dyn = 1
    elif m_n < 30:
      Z_X_dyn =  1.05 - 0.005  * m_n
    else:
      Z_X_dyn = 0.9
    Z_X_stat = 1
    print(f"Z_X_stat = {Z_X_stat}")
    print(f"Z_X_dyn = {Z_X_dyn}")

    Z_LVR = 1
    print(f"Z_L*Z_V*Z_R = {Z_LVR}")

    Z_W = 1
    print(f"Z_W = {Z_W}")

    sigma_HGstat = sigma_Hlim * Z_NTstat * Z_LVR * Z_W * Z_X_stat
    sigma_HGdyn = sigma_Hlim * Z_NTdyn *  Z_LVR * Z_W * Z_X_dyn
    print(f"sigma_HGstat = {sigma_HGstat}")
    print(f"sigma_HGdyn = {sigma_HGdyn}")

    S_Hstat = sigma_HGstat / sigma_Hstat
    S_Hdyn = sigma_HGdyn / sigma_Hdyn
    print(f"S_Hstat = {S_Hstat} soll > {S_Hstat_min}")
    print(f"S_Hdyn = {S_Hdyn} soll in {S_Hdyn_interval} liegen")
    if not (S_Hstat > S_Hstat_min):
      print("\033[91mstatische Gruebchensicherheit ist nicht erfuellt\033[0m")
      return False
    if not (S_Hdyn_interval[0] <= S_Hdyn <= S_Hdyn_interval[1]):
      print("\033[91mdynamische Gruebchensicherheit ist nicht erfuellt\033[0m")
      return False

    # Zahnfusstragfaehigkeit

    Y_epsilon = 0.25 + 0.75 / epsilon_alpha * m.cos(m.radians(beta))**2
    print(f"Y_epsilon = {Y_epsilon}")

    temp_epsilon_beta = epsilon_beta
    if epsilon_beta > 1:
      temp_epsilon_beta = 1

    temp_beta = beta
    if beta > 30:
      temp_beta = 30

    Y_beta = 1 - temp_epsilon_beta * temp_beta / 120
    print(f"Y_beta = {Y_beta}")

    x = x_1 if ist_ritzel else x_2
    G = rho_fP / m_n - h_fP / m_n + x # D.5.02
    print(f"G = {G}")

    z_n = (z_1 if ist_ritzel else z_2) / 0.9869 / m.cos(m.radians(beta))
    print(f"z_n = {z_n}")

    d = d_1 if ist_ritzel else d_2
    d_n = m_n * z_n # D.5.06

    d_bn = d_n * m.cos(m.radians(alpha_n))  # D.5.07

    d_a = d_a1 if ist_ritzel else d_a2
    d_an = d_n + d_a - d  # D.5.08

    alpha_an = m.degrees(m.acos(d_bn / d_an)) # degrees, D.5.09

    y_a = 1 / z_n * (m.pi / 2 + 2 * x * m.tan(m.radians(alpha_n))) + involute(alpha_n) - involute(alpha_an) # radiant, D.5.10

    alpha_Fan = alpha_an - m.degrees(y_a) # degrees, D.5.11
    print(f"alpha_Fan = {alpha_Fan}")

    E = m.pi / 4 * m_n - h_fP * m.tan(m.radians(alpha_n)) + s_pr / m.cos(m.radians(alpha_n)) - (1 - m.sin(m.radians(alpha_n))) * rho_fP / m.cos(m.radians(alpha_n)) # Glg D.5.01
    print(f"E = {E}")

    H = 2 / z_n * (m.pi / 2 - E / m_n) - m.pi / 3
    print(f"H = {H}")

    delta = m.pi / 6  # in radiant
    for i in range(5):
      delta = 2 * G / z_n * m.tan(delta) - H  # D.5.04
    print(f"delta = {delta}")

    rho_F = rho_fP + m_n * 2 * G**2 / m.cos(delta) / (z_n * m.cos(delta)**2 - 2 * G)  # D.5.13

    h_Fa = m_n * (0.5 * z_n * (m.cos(m.radians(alpha_n)) / m.cos(m.radians(alpha_Fan)) - m.cos(m.pi / 3 - delta)) + 0.5 * (rho_fP / m_n - G / m.cos(delta)))  # D.5.12
    print(f"h_Fa = {h_Fa}")

    s_Fn = m_n * (z_n * m.sin(m.pi / 3 - delta) + m.sqrt(3) * (G / m.cos(delta) - rho_fP / m_n))  # D.5.05
    print(f"s_Fn = {s_Fn}")

    Y_Fa = 6 * h_Fa / m_n * m.cos(m.radians(alpha_Fan)) / (s_Fn / m_n)**2 / m.cos(m.radians(alpha_n)) # Glg D.3.01
    print(f"Y_Fa = {Y_Fa}")

    L_a = s_Fn / h_Fa # Glg D.4.02

    q_s = s_Fn / 2 / rho_F  # Glg D.4.03
    print(f"q_s = {q_s}")

    assert(1 <= q_s < 8)  # Abschnitt D.4
    Y_Sa = (1.2 + 0.13 * L_a) * q_s ** (1 / (1.21 + 2.3 / L_a)) # Glg D.4.01
    print(f"Y_Sa = {Y_Sa}")

    Y_FS = Y_Sa * Y_Fa
    print(f"Y_FS = {Y_FS}")

    sigma_F0 = F_t / b / m_n * Y_FS * Y_epsilon * Y_beta
    print(f"sigma_F0 = {sigma_F0}")

    sigma_Fstat = sigma_F0 * K_S * K_V * K_Fbeta * K_Falpha
    sigma_Fdyn = sigma_F0 * K_A * K_V * K_Fbeta * K_Falpha
    print(f"sigma_Fstat = {sigma_Fstat}")
    print(f"sigma_Fdyn = {sigma_Fdyn}")

    # zulaessige Beanspruchung

    Y_NT_stat = 2.5
    Y_NT_dyn = 1
    print(f"Y_NT_stat = {Y_NT_stat}")
    print(f"Y_NT_dyn = {Y_NT_dyn}")

    Y_S = Y_Sa * (0.6 + 0.4 * epsilon_alpha)
    print(f"Y_S = {Y_S}")

    Y_delta_rel_T_stat = 0.44 * Y_S + 0.12  # Glg 5.15
    if q_s >= 1.5:
      Y_delta_rel_T_dyn = 1
    else:
      Y_delta_rel_T_dyn = 0.95
    print(f"Y_delta_rel_T_stat = {Y_delta_rel_T_stat}")
    print(f"Y_delta_rel_T_dyn = {Y_delta_rel_T_dyn}")

    # Tabelle 5.1
    Y_X_stat = 1
    if m_n <= 5:
      Y_X_dyn = 1
    elif m_n < 25:
      Y_X_dyn = 1.05 - 0.01 * m_n
    else:
      Y_X_dyn = 0.8
    print(f"Y_X_stat = {Y_X_stat}")
    print(f"Y_X_dyn = {Y_X_dyn}")

    Y_R_rel_T_stat = 1
    Y_R_rel_T_dyn = 1
    print(f"Y_R_rel_T_stat = {Y_R_rel_T_stat}")
    print(f"Y_R_rel_T_dyn = {Y_R_rel_T_dyn}")

    sigma_FGstat = sigma_FE * Y_NT_stat * Y_delta_rel_T_stat * Y_R_rel_T_stat * Y_X_stat
    sigma_FGdyn = sigma_FE * Y_NT_dyn * Y_delta_rel_T_dyn * Y_R_rel_T_dyn * Y_X_dyn
    print(f"sigma_FGstat = {sigma_FGstat}")
    print(f"sigma_FGdyn = {sigma_FGdyn}")

    print(f"S_Hstat = {S_Hstat} soll > {S_Hstat_min}")
    print(f"S_Hdyn = {S_Hdyn} soll in {S_Hdyn_interval} liegen")

    S_Fstat = sigma_FGstat / sigma_Fstat
    S_Fdyn = sigma_FGdyn / sigma_Fdyn
    print(f"S_Fstat = {S_Fstat} soll > {S_Fstat_min}")
    print(f"S_Fdyn = {S_Fdyn} soll in {S_Fdyn_interval} liegen")

    if not (S_Hstat > S_Hstat_min):
      print("\033[91mstatische Gruebchensicherheit ist nicht erfuellt\033[0m")
      return False
    if not (S_Hdyn_interval[0] <= S_Hdyn <= S_Hdyn_interval[1]):
      print("\033[91mdynamische Gruebchensicherheit ist nicht erfuellt\033[0m")
      return False

    if not (S_Fstat > S_Fstat_min):
      print("\033[91mstatische Zahnbruchsicherheit ist nicht erfuellt\033[0m")
      return False
    if not (S_Fdyn_interval[0] <= S_Fdyn <= S_Fdyn_interval[1]):
      print("\033[91mdynamische Zahnbruchsicherheit ist nicht erfuellt\033[0m")
      return False

    print()
    return True

  print()
  print("Festigkeitswerte Ritzel")
  ritzel = Festigkeitsnachweis(verzahnungs_qualitaten[0], 0, 0, True)
  print("Festigkeitswerte Rad")
  rad = Festigkeitsnachweis(verzahnungs_qualitaten[1], 0, 0, False)

  return ritzel and rad

class Getriebe:
    """
    Parameter:
    - verzahnungsqualitäten: Verzahnungsqualitäten nach DIN 3961

    Zusätzliche Parameter:
    - b: Breite
    - b_d_1_verhältnis: b/d_1 Verhältnis
    """
    def __init__(self,
            m_n : float,
            z: tuple[int, int],
            x: tuple[float, float],
            bezugsprofil : Profil,
            beta : float,
            k : int,
            **kwargs):

        self.m_n = m_n
        self.z = z
        self.x = x
        self.beta = beta
        self.k = k

        print("Getriebegeometrie")
        [print(f"{key} = {value}") for key, value in vars(self).items()]

        self.alpha_n = bezugsprofil.alpha_n
        self.h_aP = bezugsprofil.h_aP_s * self.m_n
        self.h_fP = bezugsprofil.h_fP_s * self.m_n
        self.rho_fP = bezugsprofil.rho_fP_s * self.m_n
        print(f"alpha_n = {self.alpha_n}")
        print(f"h_aP = {self.h_aP}")
        print(f"h_fP = {self.h_fP}")
        print(f"rho_fP = {self.rho_fP}")

        self.u = self.z[Rad] / self.z[Ritzel]
        print(f"u = {self.u}")

        def d(idx):
            return self.z[idx] * self.m_n / m.cos(m.radians(self.beta))
        self.d = d(Ritzel), d(Rad)
        print(f"d = {self.d}")

        if "b" in kwargs:
            self.b : float = kwargs.pop("b")
        elif "b_d_1_verhältnis" in kwargs:
            self.b : float = self.d[Ritzel] * kwargs.pop("b_d_1_verhältnis")
        else:
            raise ValueError("b or b_d_1_verhältnis must be specified as argument")
        print(f"b = {self.b}")

        self.m_t = m_n / m.cos(m.radians(self.beta))
        print(f"m_t = {self.m_t}")

        self.alpha_t = m.degrees(m.atan(m.tan(m.radians(self.alpha_n)) / m.cos(m.radians(self.beta))))
        print(f"alpha_t = {self.alpha_t}")

        def d_b(idx):
            return self.d[idx] * m.cos(m.radians(self.alpha_t))
        self.d_b = d_b(Ritzel), d_b(Rad)
        print(f"d_b = {self.d_b}")

        def d_a(idx):
            return self.d[idx] + 2 * (self.x[idx] * self.m_n + self.h_aP + self.k * self.m_n)
        self.d_a = d_a(Ritzel), d_a(Rad)
        print(f"d_a = {self.d_a}")

        def d_f(idx):
            return self.d[idx] - 2 * (self.h_fP - self.x[idx] * self.m_n)
        self.d_f = d_f(Ritzel), d_f(Rad)
        print(f"d_f = {self.d_f}")

        self.h = self.h_aP + k * m_n + self.h_fP
        print(f"h = {self.h}")

        self.alpha_wt = inverse_involute( involute(self.alpha_t) + 2 * sum(self.x) / sum(self.z) * m.tan(m.radians(self.alpha_n)) )
        print(f"alpha_wt = {self.alpha_wt}")

        def d_w(idx):
            return self.d_b[idx] / m.cos(m.radians(self.alpha_wt))
        self.d_w = d_w(Ritzel), d_w(Rad)
        print(f"d_w1 = {self.d_w}")

        self.a_w = sum(self.d_w) / 2
        print(f"a_w = {self.a_w}")

        # Profilueberdeckung

        self.p_n = self.m_n * m.pi
        print(f"p_n = {self.p_n}")

        self.p_t = self.p_n / m.cos(m.radians(self.beta))
        print(f"p_t = {self.p_t}")

        self.epsilon_alpha = (m.sqrt(self.d_a[Ritzel]**2 - self.d_b[Ritzel]**2) + m.sqrt(self.d_a[Rad]**2 - self.d_b[Rad]**2) - sum(self.d_b) * m.tan(m.radians(self.alpha_wt))) / (2 * self.p_t * m.cos(m.radians(self.alpha_t)))
        print(f"epsilon_alpha = {self.epsilon_alpha}")

        self.epsilon_beta = self.b * m.sin(m.radians(self.beta)) / self.p_n
        print(f"epsilon_beta = {self.epsilon_beta}")

        self.epsilon_gamma = self.epsilon_alpha + self.epsilon_beta
        print(f"epsilon_gamma = {self.epsilon_gamma}")
        assert self.epsilon_gamma > 1

        # Unterschnitt

        self.h_aP0 = self.h_fP
        print(f"h_aP0 = {self.h_aP0}")

        def z_min(idx):
            return 2 * m.cos(m.radians(self.beta)) * (self.h_aP0 / self.m_n - self.x[idx]) / m.sin(m.radians(self.alpha_t))**2
        self.z_min = z_min(Ritzel), z_min(Rad)
        print(f"z_min = {self.z_min}")
        assert self.z[Ritzel] > self.z_min[Ritzel]
        assert self.z[Rad] > self.z_min[Rad]

        # Spitzwerden

        def gamma(idx):
            return inverse_involute( m.pi / 2 / self.z[idx] + 2 * self.x[idx] / self.z[idx] * m.tan(m.radians(self.alpha_n)) + involute(self.alpha_t) )
        self.gamma = gamma(Ritzel), gamma(Rad)
        print(f"gamma = {self.gamma}")

        def d_amax(idx):
            return self.m_t * self.z[idx] * m.cos(m.radians(self.alpha_t)) / m.cos(m.radians(self.gamma[idx]))
        self.d_amax = d_amax(Ritzel), d_amax(Rad)
        print(f"d_amax = {self.d_amax}")
        assert self.d_a[Ritzel] <= self.d_amax[Ritzel]
        assert self.d_a[Rad] <= self.d_amax[Rad]


        assert len(kwargs) == 0, f"unknown arguments: {kwargs.keys()}"
        print()
        return
    
    def Zahnrad(self,
                idx : int,
                P : float,
                n : float,
                verzahnungsqualitäten : int | tuple[int, int],
                werkstoffe : Werkstoff | tuple[Werkstoff, Werkstoff],
                K_A : float,
                K_S : float,
                **kwargs):
        """Für zusätzliche Parameter siehe DIN_3990_11"""

        vq = (verzahnungsqualitäten, verzahnungsqualitäten) if isinstance(verzahnungsqualitäten, int) else verzahnungsqualitäten
        ws = (werkstoffe, werkstoffe) if isinstance(werkstoffe, Werkstoff) else werkstoffe
        kwargs.update((key, (value[idx] if isinstance(value, tuple) else value)) for key, value in vars(self).items())
        return DIN_3990_11(self, idx, P = P, n = n, verzahnungsqualitäten = vq, werkstoffe = ws, K_A = K_A, K_S = K_S, **kwargs)

    def Zahnräder(self,
                P : float,
                n_1 : float,
                verzahnungsqualitäten : int | tuple[int, int],
                werkstoffe : Werkstoff | tuple[Werkstoff, Werkstoff],
                K_A : float | tuple[float, float],
                K_S : float | tuple[float, float],
                **kwargs):
        """Für zusätzliche Parameter siehe DIN_3990_11"""
        
        K_A = (K_A, K_A) if isinstance(K_A, float) else K_A
        K_S = (K_S, K_S) if isinstance(K_S, float) else K_S
        kwargs1 = {key : (value[Ritzel] if isinstance(value, tuple) else value) for key, value in kwargs.items()}
        kwargs2 = {key : (value[Rad] if isinstance(value, tuple) else value) for key, value in kwargs.items()}
        return self.Zahnrad(Ritzel, P, n_1, verzahnungsqualitäten, werkstoffe, K_A[Ritzel], K_S[Ritzel], **kwargs1), self.Zahnrad(Rad, P, n_1/self.u, vq[Rad], ws[Rad], K_A[Rad], K_S[Rad], **kwargs2)

    pass


class DIN_3990_11:
    def __init__(self,
                base : Getriebe,
                idx : int,
                P : float,
                n : float,
                verzahnungsqualitäten : tuple[int, int],
                werkstoffe : tuple[Werkstoff, Werkstoff],
                K_A : float,
                K_S : float,
                A : float = 0.023,
                bild_3_1 : str = "a",
                bild_3_2 : str = "a",
                s : float = 0,
                l : float = 0,
                stützwirkung : float = False,
                doppelschrägverzahnt : bool = False,
                f_ma : float = 0,
                **kwargs):
        """
        Parameter:
        - P: Leistung [kW]
        - n: Drehzahl [1/min]
        - verzahnungsqualitäten: siehe Tabelle 3.1
        - werkstoffe
        - K_A: siehe Tabelle A.1
        - K_S: ersetz K_A für die statische Berechnung
        - A: siehe Tabelle 3.2
        - bild_3_1: a-f, siehe Bild 3.1
        - bild_3_2: a-e, siehe Bild 3.2
        - s: siehe Bild 3.2
        - l: siehe Bild 3.2
        - stützwirkung: siehe Bild 3.2
        - d_sh1: Wellendurchmesser
        - doppelschrägverzahnt
        - f_ma:  Flankenlinien-Herstellabweichung, siehe Abschnitt 3.4.2.4
        """

        self.__dict__.update((key, value) for key, value in locals().items() if key not in ("base", "idx", "kwargs"))
        # self.P = P
        # self.n = n
        # self.verzahnungsqualität = verzahnungsqualität
        # self.werkstoff = werkstoff
        # self.K_A = K_A
        # self.K_S = K_S

        assert(all(q in range(6, 13) for q in self.verzahnungsqualität))

        print("DIN 3990-11 " + ("Ritzel" if idx == Ritzel else "Rad"))
        print()

        print("Parameter")
        [print(f"{key} = {value}") for key, value in vars(self).items()]

        self.__dict__.update(kwargs)
        self.base = base
        #self.idx = idx

        if self.doppelschrägverzahnt:
            self.d_B = self.b / 2
            print(f"d_B = {self.d_B}")

        self.v = self.n / 60 * 2 * m.pi * self.d / 2000
        print(f"v = {self.v}")

        self.F_t = 1000 * self.P / self.v
        print(f"F_t = {self.F_t}")
  
        self.T = self.F_t * self.d / 2000
        print(f"T = {self.T}")

        assert(self.F_t / self.b * self.K_A >= 100)  # Abschnitt 3.4.1
        self.F_m = self.F_t * self.K_A * self.K_V
        print(f"F_m = {self.F_m}")


        return

    def __getattr__(self, name):
        funcname = "_" + name
        assert hasattr(self, funcname)
        value = getattr(self, funcname)()
        setattr(self, name, value)
        print(f"{name} = {value}")
        return value

    def _K_V(self):
        temp1 = self.base.z[Ritzel] * self.v / 100 * m.sqrt(self.u**2 / (1 + self.u**2))
        assert(temp1 < 10)  # Glg 3.04

        temp2 = max(self.K_A * self.F_t / self.b, 100)

        def K_V(geradverzahnt):
            if geradverzahnt:
                match self.verzahnungsqualität:
                    case 6:
                        K_1 = 9.6
                    case 7:
                        K_1 = 15.3
                    case 8:
                        K_1 = 24.5
                    case 9:
                        K_1 = 34.5
                    case 10:
                        K_1 = 53.6
                    case 11:
                        K_1 = 76.6
                    case 12:
                        K_1 = 122.5
                K_2 = 0.0193
            else:
                match self.verzahnungsqualität:
                    case 6:
                        K_1 = 8.5
                    case 7:
                        K_1 = 13.6
                    case 8:
                        K_1 = 21.8
                    case 9:
                        K_1 = 30.7
                    case 10:
                        K_1 = 47.7
                    case 11:
                        K_1 = 68.2
                    case 12:
                        K_1 = 109.1
                K_2 = 0.0087
            return 1 + (K_1 / temp2 + K_2) * temp1

        if self.beta == 0:
            return K_V(True)
        elif self.epsilon_beta >= 1:
            return K_V(False)
        else:
            K_Valpha = K_V(True)
            K_Vbeta = K_V(False)
            return K_Valpha - self.epsilon_beta * (K_Valpha - K_Vbeta)
    
    def _F_betay(self):
        """Glg 0.08"""
        return self.F_betax - self.gamma_beta
    def _F_betax(self):
        """Glg 0.09"""
        match self.bild_3_1:
            case "a", "f":
                multi = -1
            case "b", "e":
                multi = 1
            case "c":
                B_s = 1.5 if self.doppelschrägverzahnt else 1
                multi = 1 if abs(self.K_s) * self.l * self.s / self.base.d[Ritzel]**2 * (self.base.d[Ritzel] / self.d_sh1)**4 <= B_s else -1
            case "d":
                B_s = 1.5 if self.doppelschrägverzahnt else 1
                multi = 1 if abs(self.K_s) * self.l * self.s / self.base.d[Ritzel]**2 * (self.base.d[Ritzel] / self.d_sh1)**4 >= B_s - 0.3 else -1
            case _:
                raise ValueError

        F_betax = abs(1.33 * self.f_sh + multi * self.f_ma)
        assert F_betax >= self.F_betaxmin
        return F_betax
    def _F_betaxmin(self):
        return max(0.005 * self.F_m / self.b, 0.5 * self.f_Hbeta)
    def _f_sh(self):
        if self.s == 0:
            temp = 1
        else:
            temp = self.K_s - self.l * self.s / self.base.d[Ritzel]**2 * (self.base.d[Ritzel] / self.d_sh)**4

        if not self.doppelschrägverzahnt:   # Glg 3.14
            return self.F_m / self.b * self.A * (abs(1 + temp - 0.3) + 0.3) * (self.b / self.base.d[Ritzel])**2
        else:   # Glg 3.15
            return self.F_m / self.b * 2 * self.A * (abs(1.5 + temp - 0.3) + 0.3) * (self.b_B / self.base.d[Ritzel])**2
    def _y_beta(self):
        def y_beta(idx):
            if self.werkstoffe[idx].art == Werkstoff.Art.EH_Stahl:
                return 
        y_beta1 = y_beta(Ritzel)
        y_beta2 = y_beta(Rad)
        return (y_beta1 + y_beta2) / 2

    pass

"""

berechnung von Ritzel und Rad kann nicht getrennt durchgeführt werden

"""

getriebe = Getriebe(m_n = 4,
                    z = (19, 104),
                    x = (0.5, 0.15),
                    bezugsprofil = Normalprofil1,
                    beta = 0,
                    k = 0,
                    b_d_1_verhältnis = 0.64)

werkstoff = Werkstoff(Werkstoff.Art.EH_Stahl, 1500, 860)
ritzel, rad = getriebe.Zahnräder(P = 55,
            n = 980,
            verzahnungsqualitäten = (6, 7),
            werkstoffe = (werkstoff, werkstoff),
            K_A = 1.75,
            K_S = 2.5)


DIN_3990_11()

assert(self.b / self.m_n <= 30) # Konstruktionsvorgaben, Tabelle 4

sys.exit()

execute(P = 55,
            n = 980,
            b_d_1_verhaeltnis = 0.64,
            verzahnungs_qualitaten = (6, 7),
            m_n = 4,
            z_1 = 19,
            z_2 = 104,
            x_1 = 0.5,
            x_2 = 0.15,
            alpha_n = 20,
            beta = 0,
            k = 0,
            sigma_Hlim = 1500,
            sigma_FE = 860,
            K_A = 1.75,
            K_S = 2.5)


"""K_A_in = (1.5, 1.6, 1.75, 1.85)
m_n_in = (3, 4, 5, 6, 7)
K_S_in = range(20, 31)
_print = print
print = lambda *args, **kwargs: None
x_2 = 0.15
i = 0
hits = 0
#for K_S in K_S_in:
for K_S in K_S_in:
  K_S = K_S / 10
  for z_1 in range(19, 20):#range(14, 20):
    x_1min = (1.25 - z_1 * m.sin(m.radians(alpha_n))**2/2)*100
    for z_2 in range(104, 105):#range(5 * z_1, m.ceil(5.5 * z_1)):
      if (z_2 / z_1).is_integer():
        continue
      for x_1 in range(m.ceil(x_1min), 150):
        x_1 = x_1 / 100
        for K_A in K_A_in:
          for m_n in m_n_in:
            i += 1

            try:
              res = func(P = 55,
                  n = 980,
                  b_d_1_verhaeltnis = 0.64,
                  m_n = m_n,
                  z_1 = z_1,
                  z_2 = z_2,
                  x_1 = x_1,
                  x_2 = x_2,
                  alpha_n = 20,
                  beta = 0,
                  k = 0,
                  sigma_Hlim = 1500,
                  sigma_FE = 860,
                  K_A = K_A,
                  K_S = K_S)
              if res:
                hits += 1
                _print(f"i = {i}, K_A = {K_A}, K_S = {K_S}, m_n = {m_n}, z_1 = {z_1}, z_2 = {z_2}, x_1 = {x_1}, x_2 = {x_2}")
                #sys.exit()
            except AssertionError:
              pass
        
_print(f"iterations: ~{i}")
_print(f"hits: {hits}")"""
