from dataclasses import dataclass
from enum import Enum, StrEnum
import math as m, sys, multiprocessing, copy
from scipy import optimize

S_Hstat_min = 1.3
S_Hdyn_interval = (1.2, 1.5)
S_Fstat_min = 3.5
S_Fdyn_interval = (1.5, 2)


class HiddenPrints:
    def __init__(self, switch : bool):
        self.switch = switch

    def __enter__(self):
        if self.switch:
            self._original_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.switch:
            sys.stdout.close()
            sys.stdout = self._original_stdout

def to_float(val) -> float:
    return val

def involute(alpha):
    return m.tan(m.radians(alpha)) - m.radians(alpha)
def inverse_involute(alpha, anfangswert = 20):
    try:
        return float(optimize.newton(lambda x: involute(x) - alpha, anfangswert))
    except RuntimeError:
        assert(False)
   
Ritzel = 0
Rad = 1

@dataclass
class Profil:
    alpha_n : float
    h_aP_s : float
    h_fP_s : float
    rho_fP_s : float
Normalprofil1 =     Profil(20, 1, 1.25, 0.250)
Normalprofil2 =     Profil(20, 1, 1.25, 0.375)
Protuberanzprofil = Profil(20, 1, 1.40, 0.400)

@dataclass
class Werkstoff:
    class Art(Enum):
        Baustahl = 0
        vergüteterStahl = 1
        einsatzgehärteterStahl = 2
        randschichtgehärteterStahl = 3
        nitrierterStahl = 4
        nitrokarburierterStahl = 5
        
    art : Art
    sigma_Hlim : float
    sigma_FE : float
    HB : float
    """Nur das Intervall (130; 470) ist relevant"""

class Tabelle_3_2(Enum):
    ohneBreitenballigkeltOderEndrücknahme = 0
    mitSinnvollerBreitenballigkeit = 1
    mitSinnvollerEndrücknahme = 2

class Bild_3_1(Enum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4
    f = 5

class Bild_3_2(Enum):
    a = 0
    b = 1
    c = 2
    d = 3
    e = 4

class Fertigungsverfahren(Enum):
    wälzgefrästWälzgestoßenWälzgehobelt = 0
    geläpptGeschliffenGeschabt = 1

class Geometrie:
    def __init__(self,
            m_n : float,
            z: tuple[int, int],
            x: tuple[float, float],
            bezugsprofil : Profil,
            beta : float,
            k : int,
            **kwargs):
        """
        Parameter:

        Zusätzliche Parameter:
        - b: Breite
        - b_d_1_verhältnis: b/d_1 Verhältnis
        """

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
            self.b = to_float(kwargs.pop("b"))
        elif "b_d_1_verhältnis" in kwargs:
            self.b = self.d[Ritzel] * to_float(kwargs.pop("b_d_1_verhältnis"))
        else:
            raise ValueError("b or b_d_1_verhältnis must be specified as argument")
        print(f"b = {self.b}")

        self.m_t = m_n / m.cos(m.radians(self.beta))
        print(f"m_t = {self.m_t}")

        self.alpha_t = m.degrees(m.atan(m.tan(m.radians(self.alpha_n)) / m.cos(m.radians(self.beta))))
        print(f"alpha_t = {self.alpha_t}")

        self.beta_b = m.degrees(m.atan(m.tan(m.radians(self.beta)) * m.cos(m.radians(self.alpha_t))))
        print(f"beta_b = {self.beta_b}")

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
  
class Getriebe:
    def __init__(self,
                geometrie : Geometrie,
                P : float,
                n_1 : float,
                verzahnungsqualität : tuple[int, int],
                werkstoff : tuple[Werkstoff, Werkstoff],
                K_A : float,
                K_S : float,
                R_z: tuple[float, float],
                doppelschrägverzahnt: bool = False,
                innenverzahnt: bool = False,
                s_pr: float = 0,
                output = True,
                **kwargs):
        """
        Parameter:
        - geometrie
        - P: Leistung [kW]
        - n_1: Antriebsdrehzahl [1/min]
        - verzahnungsqualität: siehe Tabelle 3.1
        - werkstoff
        - K_A: siehe Tabelle A.1
        - K_S: ersetz K_A für die statische Berechnung
        - R_z: gemittelte Rauhtiefe
        - doppelschrägverzahnt
        - innenverzahnt
        - s_pr: float, Fussfreischnitt

        Zusätzliche Parameter:
        - K_Hbeta : tuple[float, float]
          or 
            - A: Tabelle_3_2, siehe Tabelle 3.2
            - f_ma: tuple[float, float], siehe Abschnitt 3.4.2.4
              if f_ma != 0:
                - bild_3_1: tuple[Bild_3_1, Bild_3_1], siehe Bild 3.1
            - s: tuple[float, float], siehe Bild 3.2
              if s != 0:
                - stützwirkung: tuple[bool, bool], siehe Bild 3.2
                - bild_3_2: tuple[Bild_3_2, Bild_3_2], siehe Bild 3.2
                - l: tuple[float, float], siehe Bild 3.2
                - d_sh: tuple[float, float], Wellendurchmesser
        - Z_LVRdyn: float, siehe Abschnitt 4.8
          or
            - fertigungsverfahren: tuple[Fertigungsverfahren, Fertigungsverfahren], siehe Abschnitt 4.8
        - alle Werte die von einem if not hasattr umgeben sind
        """
        
        with HiddenPrints(not output):
            print("Parameter")
            [print(f"{key} = {value}") for key, value in locals().items() if key not in ("self", "kwargs", "output")]
            print()

            print("Zusätzliche Parameter")
            [print(f"{key} = {value}") for key, value in kwargs.items()]
            print()

            def to_tuple(val):
                return val if isinstance(val, tuple) else (val, val)

            self.geometrie = geometrie
            self.P = P
            self.n = (n_1, n_1 / self.geometrie.u)
            self.verzahnungsqualität = to_tuple(verzahnungsqualität)
            self.werkstoff = to_tuple(werkstoff)
            self.K_A = K_A
            self.K_S = K_S
            self.R_z = R_z
            self.doppelschrägverzahnt = doppelschrägverzahnt
            self.innenverzahnt = innenverzahnt
            self.s_pr = s_pr
            self.__dict__.update(kwargs)

            self._init()
        return

    def _K_V(self):
        # Glg 3.04
        temp1 = self.geometrie.z[Ritzel] * self.v / 100 * m.sqrt(self.geometrie.u**2 / (1 + self.geometrie.u**2))
        assert temp1 < 10  

        temp2 = max(self.K_A * self.F_t / self.geometrie.b, 100)

        def K_V(idx, geradverzahnt) -> float:
            if geradverzahnt:
                match self.verzahnungsqualität[idx]:
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
                match self.verzahnungsqualität[idx]:
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
        if self.geometrie.beta == 0:
            return K_V(Ritzel, True), K_V(Rad, True)
        elif self.geometrie.epsilon_beta >= 1:
            return K_V(Ritzel, False), K_V(Rad, False)
        else:
            K_Valpha = K_V(Ritzel, True), K_V(Rad, True)
            K_Vbeta = K_V(Ritzel, False), K_V(Rad, False)
            print(f"K_Valpha = {K_Valpha}")
            print(f"K_Vbeta = {K_Vbeta}")
            def interp(idx) -> float:
                return K_Valpha[idx] - self.geometrie.epsilon_beta * (K_Valpha[idx] - K_Vbeta[idx])
            return interp(Ritzel), interp(Rad)

    def _K_Hbeta(self):
        # Abschnitt 3.4.1
        assert(self.F_t / self.geometrie.b * self.K_A >= 100)

        # Glg 3.07
        def F_m(idx) -> float:
            return self.F_t * self.K_A * self.K_V[idx]
        self.F_m = F_m(Ritzel), F_m(Rad)
        print(f"F_m = {self.F_m}")
        
        def K_s(idx):
            stütz = self.stützwirkung[idx]
            match self.bild_3_2[idx]:
                case Bild_3_2.a:
                    return 0.48 if stütz else 0.8
                case Bild_3_2.b:
                    return -0.48 if stütz else -0.8
                case Bild_3_2.c:
                    return 1.33
                case Bild_3_2.d:
                    return -0.36 if stütz else -0.6
                case Bild_3_2.e:
                    return -0.6 if stütz else -1.0
            raise ValueError(f"Unknown bild_3_2 option {self.bild_3_2[idx]}")
        
        # Glg 3.14, 3.15
        def f_sh(idx) -> float:
            def temp(idx) -> float:
                return 0. if self.s[idx] == 0 else K_s(idx) * self.l[idx] * self.s[idx] / self.geometrie.d[Ritzel]**2 * (self.geometrie.d[Ritzel] / self.d_sh[idx])**4
            if self.A == Tabelle_3_2.ohneBreitenballigkeltOderEndrücknahme:
                A = 0.023
            elif self.A == Tabelle_3_2.mitSinnvollerBreitenballigkeit:
                A = 0.012
            elif self.A == Tabelle_3_2.mitSinnvollerEndrücknahme:
                A = 0.016

            if not self.doppelschrägverzahnt:
                return self.F_m[idx] / self.geometrie.b * A * (abs(1 + temp(idx) - 0.3) + 0.3) * (self.geometrie.b / self.geometrie.d[Ritzel])**2
            else:
                return self.F_m[idx] / self.geometrie.b * 2 * A * (abs(1.5 + temp(idx) - 0.3) + 0.3) * (self.b_B / self.geometrie.d[Ritzel])**2
        self.f_sh = f_sh(Ritzel), f_sh(Rad)
        print(f"f_sh = {self.f_sh}")

        # Glg 3.09
        def F_betax(idx) -> float:
            def multi(idx):
                match self.bild_3_1[idx]:
                    case Bild_3_1.a | Bild_3_1.f:
                        return -1
                    case Bild_3_1.b | Bild_3_1.e:
                        return 1
                    case Bild_3_1.c:
                        B_s = 1.5 if self.doppelschrägverzahnt else 1
                        return 1 if abs(self.K_s) * self.l * self.s / self.d[Ritzel]**2 * (self.d[Ritzel] / self.d_sh[Ritzel])**4 <= B_s else -1
                    case Bild_3_1.d:
                        B_s = 1.5 if self.doppelschrägverzahnt else 1
                        return 1 if abs(self.K_s) * self.l * self.s / self.d[Ritzel]**2 * (self.d[Ritzel] / self.d_sh[Ritzel])**4 >= B_s - 0.3 else -1
                    case _:
                        raise ValueError(f"Unknown bild_3_1 option {self.bild_3_1[idx]}")
            if self.f_ma[idx] == 0:
                return abs(1.33 * self.f_sh[idx])
            else:
                return abs(1.33 * self.f_sh[idx] + multi(idx) * self.f_ma[idx])
        self.F_betax = F_betax(Ritzel), F_betax(Rad)
        print(f"F_betax = {self.F_betax}")

        # Abschnitt 3.4.2.6
        def y_beta(idx : int):
            werkstoff = self.werkstoff[idx]
            match werkstoff.art:
                case Werkstoff.Art.Baustahl | Werkstoff.Art.vergüteterStahl:
                    y_beta = 320 / werkstoff.sigma_Hlim * self.F_betax[idx]
                    if self.v <= 5:
                        pass
                    elif self.v <= 10:
                        assert y_beta <= 25600 / werkstoff.sigma_Hlim
                    else:
                        assert y_beta <= 12800 / werkstoff.sigma_Hlim
                    return y_beta
                case Werkstoff.Art.einsatzgehärteterStahl | Werkstoff.Art.nitrierterStahl | Werkstoff.Art.nitrokarburierterStahl:
                    y_beta = 0.15 * self.F_betax[idx]
                    assert y_beta <= 6
                    return y_beta
            raise NotImplementedError
        self.y_beta = (y_beta(Ritzel) + y_beta(Rad)) / 2
        print(f"y_beta = {self.y_beta}")
        
        # Glg 3.08
        def F_betay(idx):
            return self.F_betax[idx] - self.y_beta
        self.F_betay = F_betay(Ritzel), F_betay(Rad)
        print(f"F_betay = {self.F_betay}")

        self.c_gamma = 20
        
        # Glg 3.20, 3.21
        def K_Hbeta(idx):
            val = 1 + self.c_gamma * self.F_betay[idx] / (2 * self.F_m[idx] / self.geometrie.b)
            if val <= 2:
                return val
            return m.sqrt(2 * self.c_gamma * self.F_betay[idx] / (self.F_m[idx] / self.geometrie.b))
        return K_Hbeta(Ritzel), K_Hbeta(Rad)

    def _K_Fbeta(self):
        # Abschnitt 3.4.1
        assert(self.F_t / self.geometrie.b * self.K_A >= 100)

        # Glg 3.22
        def K_Fbeta(idx):
            h_b = min(self.geometrie.h / self.geometrie.b, 1. / 3.)
            return m.pow(self.K_Hbeta[idx], (1 / (1 + h_b + h_b**2)))
        return K_Fbeta(Ritzel), K_Fbeta(Rad)

    def _init(self):
        assert all(vq in range(6, 13) for vq in self.verzahnungsqualität)
        
        print("DIN 3990-11")
        print(f"n = {self.n}")

        if self.doppelschrägverzahnt:
            self.b_B = self.geometrie.b / 2
            print(f"d_B = {self.d_B}")

        self.v = self.n[Ritzel] / 60 * 2 * m.pi * self.geometrie.d[Ritzel] / 2000
        print(f"v = {self.v}")

        self.F_t = 1000 * self.P / self.v
        print(f"F_t = {self.F_t}")
  
        def T(idx):
            return self.F_t * self.geometrie.d[idx] / 2000
        self.T = T(Ritzel), T(Rad)
        print(f"T = {self.T}")

        if not hasattr(self, "K_V"):
            self.K_V = self._K_V()
        print(f"K_V = {self.K_V}")

        if not hasattr(self, "K_Hbeta"):
            self.K_Hbeta = self._K_Hbeta()
        print(f"K_Hbeta = {self.K_Hbeta}")
        
        if not hasattr(self, "K_Fbeta"):
            self.K_Fbeta = self._K_Fbeta()
        print(f"K_Fbeta = {self.K_Fbeta}")
        
        # Abschnitt 4.5
        if self.geometrie.beta == 0:
            self.Z_epsilon = m.sqrt((4. - self.geometrie.epsilon_alpha) / 3.)
        elif self.geometrie.epsilon_beta >= 1.:
            self.Z_epsilon = m.sqrt(1 / self.geometrie.epsilon_alpha)
        else:
            self.Z_epsilon = m.sqrt((4. - self.geometrie.epsilon_alpha) / 3. * (1. - self.geometrie.epsilon_beta) + self.geometrie.epsilon_beta / self.geometrie.epsilon_alpha)
        print(f"Z_epsilon = {self.Z_epsilon}")

        # Glg 5.09
        self.Y_epsilon = 0.25 + 0.75 / self.geometrie.epsilon_alpha * m.cos(m.radians(self.geometrie.beta_b))**2
        print(f"Y_epsilon = {self.Y_epsilon}")

        # Tabelle 3.3
        def K_Halpha_und_Falpha(idx : int):
            linienbelastung = self.F_t / self.geometrie.b * self.K_A
            art = self.werkstoff[idx].art
            qualität = self.verzahnungsqualität[idx]
            Art = Werkstoff.Art

            if art in (Art.einsatzgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl):
                if self.geometrie.beta == 0:
                    if linienbelastung > 100:
                        match qualität:
                            case 6 | 7:
                                return 1., 1.
                            case 8:
                                return 1.1, 1.1
                            case 9:
                                return 1.2, 1.2
                            case 10 | 11 | 12:
                                K_H = 1 / self.Z_epsilon**2
                                K_F = 1 / self.Y_epsilon**2
                                assert K_H >= 1.2
                                assert K_F >= 1.2
                                return K_H, K_F
                    else:
                        if qualität <= 6:
                            K_H = 1 / self.Z_epsilon**2
                            K_F = 1 / self.Y_epsilon**2
                            assert K_H >= 1.2
                            assert K_F >= 1.2
                            return K_H, K_F
                else:
                    if linienbelastung > 100:
                        match qualität:
                            case 6:
                                return 1., 1.
                            case 7:
                                return 1.1, 1.1
                            case 8:
                                return 1.2, 1.2
                            case 9:
                                return 1.4, 1.4
                            case 10 | 11 | 12:
                                K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                                assert K >= 1.4
                                return K, K
                    else:
                        if qualität <= 6:
                            K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                            assert K >= 1.4
                            return K, K
            else:
                if self.geometrie.beta == 0:
                    if linienbelastung > 100:
                        match qualität:
                            case 6 | 7 | 8:
                                return 1., 1.
                            case 9:
                                return 1.1, 1.1
                            case 10:
                                return 1.2, 1.2
                            case 11 | 12:
                                K_H = 1 / self.Z_epsilon**2
                                K_F = 1 / self.Y_epsilon**2
                                assert K_H >= 1.2
                                assert K_F >= 1.2
                                return K_H, K_F
                    else:
                        if qualität <= 6:
                            K_H = 1 / self.Z_epsilon**2
                            K_F = 1 / self.Y_epsilon**2
                            assert K_H >= 1.2
                            assert K_F >= 1.2
                            return K_H, K_F
                else:
                    if linienbelastung > 100:
                        match qualität:
                            case 6 | 7:
                                return 1., 1.
                            case 8:
                                return 1.1, 1.1
                            case 9:
                                return 1.2, 1.2
                            case 10:
                                return 1.4, 1.4
                            case 11 | 12:
                                K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                                assert K >= 1.4
                                return K, K
                    else:
                        if qualität <= 6:
                            K = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
                            assert K >= 1.4
                            return K, K

            raise ValueError
        K_ritzel = K_Halpha_und_Falpha(Ritzel)
        K_rad = K_Halpha_und_Falpha(Rad)
        if not hasattr(self, "K_Halpha"):
            self.K_Halpha = K_ritzel[0], K_rad[0]
        print(f"K_Halpha = {self.K_Halpha}")
        if not hasattr(self, "K_Falpha"):
            self.K_Falpha = K_ritzel[1], K_rad[1]
        print(f"K_Falpha = {self.K_Falpha}")

        # Grübchentragfähigkeit
        
        # Glg 4.12
        self.M_1 = m.tan(m.radians(self.geometrie.alpha_wt)) / m.sqrt(
            (m.sqrt(self.geometrie.d_a[Ritzel]**2 / self.geometrie.d_b[Ritzel]**2 - 1) - 2 * m.pi / self.geometrie.z[Ritzel]) *
            (m.sqrt(self.geometrie.d_a[Rad]**2 / self.geometrie.d_b[Rad]**2 - 1) - (self.geometrie.epsilon_alpha - 1) * 2 * m.pi / self.geometrie.z[Rad]))
        print(f"M_1 = {self.M_1}")

        # Glg 4.13
        self.M_2 = m.tan(m.radians(self.geometrie.alpha_wt)) / m.sqrt(
            (m.sqrt(self.geometrie.d_a[Rad]**2 / self.geometrie.d_b[Rad]**2 - 1) - 2 * m.pi / self.geometrie.z[Rad]) *
            (m.sqrt(self.geometrie.d_a[Ritzel]**2 / self.geometrie.d_b[Ritzel]**2 - 1) - (self.geometrie.epsilon_alpha - 1) * 2 * m.pi / self.geometrie.z[Ritzel]))
        print(f"M_2 = {self.M_2}")

        # Abschnitt 4.2
        if self.geometrie.beta == 0:
            self.Z_B = max(1., self.M_1)
            self.Z_D = max(1., self.M_2)
        elif self.geometrie.epsilon_beta >= 1:
            self.Z_B = 1.
            self.Z_D = 1.
        else:
            self.Z_B = max(1., self.M_1 - self.geometrie.epsilon_beta * (self.M_1 - 1))
            self.Z_D = max(1., self.M_2 - self.geometrie.epsilon_beta * (self.M_2 - 1))
        if self.innenverzahnt:
            self.Z_D = 1.
        print(f"Z_B = {self.Z_B}")
        print(f"Z_D = {self.Z_D}")

        # Glg. 4.14
        self.Z_H = m.sqrt(2 * m.cos(m.radians(self.geometrie.beta_b)) * m.cos(m.radians(self.geometrie.alpha_wt))
                          / m.cos(m.radians(self.geometrie.alpha_t))**2 / m.sin(m.radians(self.geometrie.alpha_wt)))
        print(f"Z_H = {self.Z_H}")

        # Tabelle 4.1
        ws1, ws2 = self.werkstoff
        Art = Werkstoff.Art
        stahl = (Art.Baustahl, Art.vergüteterStahl, Art.einsatzgehärteterStahl, Art.randschichtgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl)
        if ws1.art in stahl and ws2.art in stahl:
            self.Z_E = 189.8
        else:
            raise NotImplementedError
        print(f"Z_E = {self.Z_E}")

        # Glg 4.18
        self.Z_beta = m.sqrt(m.cos(m.radians(self.geometrie.beta)))
        print(f"Z_beta = {self.Z_beta}")

        # Glg 4.20
        self.R_z100 = sum(self.R_z) / 2 * m.pow(100 / self.geometrie.a_w, 1/3)
        print(f"R_z100 = {self.R_z100}")

        # Glg 4.19
        self.Z_LVRstat = 1.
        if not hasattr(self, "Z_LVRdyn"):
            if self.fertigungsverfahren == (Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt, Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt):
                self.Z_LVRdyn = 0.85
            else:
                if self.fertigungsverfahren == (Fertigungsverfahren.geläpptGeschliffenGeschabt, Fertigungsverfahren.geläpptGeschliffenGeschabt):
                    if self.R_z100 > 4:
                        self.Z_LVRdyn = 0.92
                    else:
                        self.Z_LVRdyn = 1.
                elif Fertigungsverfahren.wälzgefrästWälzgestoßenWälzgehobelt in self.fertigungsverfahren and Fertigungsverfahren.geläpptGeschliffenGeschabt in self.fertigungsverfahren:
                    assert self.R_z100 <= 4
                    self.Z_LVRdyn = 0.92
        print(f"Z_LVRstat = {self.Z_LVRstat}")
        print(f"Z_LVRdyn = {self.Z_LVRdyn}")

        # Glg 4.23
        wsart1 = self.werkstoff[Ritzel].art
        wsart2 = self.werkstoff[Rad].art
        Art = Werkstoff.Art
        weich = (Art.Baustahl, Art.vergüteterStahl)
        hart = (Art.einsatzgehärteterStahl, Art.randschichtgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl)
        if (wsart1 in weich or wsart2 in weich) and (wsart1 in hart or wsart2 in hart) and self.R_z <= 6:
            HB = min(self.werkstoff[Ritzel].HB, self.werkstoff[Rad].HB)
            if HB < 130:
                self.Z_W = 1.2
            elif HB > 470:
                self.Z_W = 1.
            else:
                self.Z_W = 1.2 - (HB - 130) / 1700
        else:
            self.Z_W = 1.
        print(f"Z_W = {self.Z_W}")

        # Tabelle 4.2
        def Z_Xdyn(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, Art.vergüteterStahl):
                return 1.
            elif wsart in (Art.einsatzgehärteterStahl, ):
                if self.geometrie.m_n <= 10:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.05 - 0.005  * self.geometrie.m_n
                else:
                    return 0.9
            elif wsart in (Art.nitrierterStahl, Art.nitrokarburierterStahl):
                if self.geometrie.m_n <= 7.5:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.08 - 0.011  * self.geometrie.m_n
                else:
                    return 0.75
            raise NotImplementedError
        self.Z_Xstat = 1., 1.
        self.Z_Xdyn = Z_Xdyn(Ritzel), Z_Xdyn(Rad)
        print(f"Z_Xstat = {self.Z_Xstat}")
        print(f"Z_Xdyn = {self.Z_Xdyn}")

        # Tabelle 4.3
        def Z_NT(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, ):
                return 1.6, 1.0
            elif wsart in (Art.vergüteterStahl, Art.einsatzgehärteterStahl):
                return 1.6, 1.0
            elif wsart in (Art.nitrierterStahl, ):
                return 1.3, 1.0
            elif wsart in (Art.nitrokarburierterStahl, ):
                return 1.1, 1.0
            raise NotImplementedError
        Z_NTritzel = Z_NT(Ritzel)
        Z_NTrad = Z_NT(Rad)
        self.Z_NTstat = Z_NTritzel[0], Z_NTrad[0]
        self.Z_NTdyn = Z_NTritzel[1], Z_NTrad[1]
        print(f"Z_NTstat = {self.Z_NTstat}")
        print(f"Z_NTdyn = {self.Z_NTdyn}")

        # Glg 4.03
        def sigma_HGstat(idx):
            return self.werkstoff[idx].sigma_Hlim * self.Z_NTstat[idx] * self.Z_LVRstat * self.Z_W * self.Z_Xstat[idx]
        def sigma_HGdyn(idx):
            return self.werkstoff[idx].sigma_Hlim * self.Z_NTdyn[idx] * self.Z_LVRdyn * self.Z_W * self.Z_Xdyn[idx]
        self.sigma_HGstat = sigma_HGstat(Ritzel), sigma_HGstat(Rad)
        self.sigma_HGdyn = sigma_HGdyn(Ritzel), sigma_HGdyn(Rad)
        print(f"sigma_HGstat = {self.sigma_HGstat}")
        print(f"sigma_HGdyn = {self.sigma_HGdyn}")

        # Glg 4.02
        self.sigma_H0 = self.Z_H * self.Z_E * self.Z_epsilon * self.Z_beta * m.sqrt(self.F_t / self.geometrie.d[Ritzel] / self.geometrie.b * (self.geometrie.u + 1) / self.geometrie.u)
        print(f"sigma_H0 = {self.sigma_H0}")

        # Glg 4.01
        def sigma_Hstat(idx):
            return (self.Z_B if idx == Ritzel else self.Z_D) * self.sigma_H0 * m.sqrt(self.K_S * self.K_V[idx] * self.K_Hbeta[idx] * self.K_Halpha[idx])
        def sigma_Hdyn(idx):
            return (self.Z_B if idx == Ritzel else self.Z_D) * self.sigma_H0 * m.sqrt(self.K_A * self.K_V[idx] * self.K_Hbeta[idx] * self.K_Halpha[idx])
        self.sigma_Hstat = sigma_Hstat(Ritzel), sigma_Hstat(Rad)
        self.sigma_Hdyn = sigma_Hdyn(Ritzel), sigma_Hdyn(Rad)
        print(f"sigma_Hstat = {self.sigma_Hstat}")
        print(f"sigma_Hdyn = {self.sigma_Hdyn}")

        # Glg 4.11
        def S_Hstat(idx):
            return self.sigma_HGstat[idx] / self.sigma_Hstat[idx]
        def S_Hdyn(idx):
            return self.sigma_HGdyn[idx] / self.sigma_Hdyn[idx]
        self.S_Hstat = S_Hstat(Ritzel), S_Hstat(Rad)
        self.S_Hdyn = S_Hdyn(Ritzel), S_Hdyn(Rad)
        print(f"S_Hstat = {self.S_Hstat}")
        print(f"S_Hdyn = {self.S_Hdyn}")

        # Zahnfußtragfähigkeit

        # Glg D.1.01
        def z_n(idx):
            return self.geometrie.z[idx] / m.cos(m.radians(self.geometrie.beta_b))**2 / m.cos(m.radians(self.geometrie.beta))
        self.z_n = z_n(Ritzel), z_n(Rad)
        print(f"z_n = {self.z_n}")

        if not self.innenverzahnt:
            # Glg D.5.01
            self.E = (m.pi / 4 * self.geometrie.m_n - self.geometrie.h_fP * m.tan(m.radians(self.geometrie.alpha_n)) + self.s_pr / m.cos(m.radians(self.geometrie.alpha_n))
                      - (1 - m.sin(m.radians(self.geometrie.alpha_n))) * self.geometrie.rho_fP / m.cos(m.radians(self.geometrie.alpha_n))) 
            print(f"E = {self.E}")

            # Glg D.5.02
            def G(idx):
                return self.geometrie.rho_fP / self.geometrie.m_n - self.geometrie.h_fP / self.geometrie.m_n + self.geometrie.x[idx]
            self.G = G(Ritzel), G(Rad)
            print(f"G = {self.G}")

            # Glg D.5.03
            def H(idx):
                return 2 / self.z_n[idx] * (m.pi / 2 - self.E / self.geometrie.m_n) - m.pi / 3
            self.H = H(Ritzel), H(Rad)
            print(f"H = {self.H}")

            # Glg D.5.04
            def delta(idx):
                delta = m.degrees(m.pi / 6)
                for i in range(5):
                    delta = m.degrees(2 * self.G[idx] / self.z_n[idx] * m.tan(m.radians(delta)) - self.H[idx])
                return delta
            self.delta = delta(Ritzel), delta(Rad)
            print(f"delta = {self.delta}")

            # Glg D.5.05
            def s_Fn(idx):
                return self.geometrie.m_n * (self.z_n[idx] * m.sin(m.pi / 3 - m.radians(self.delta[idx])) + m.sqrt(3) * (self.G[idx] / m.cos(m.radians(self.delta[idx]))
                                                                                                                         - self.geometrie.rho_fP / self.geometrie.m_n))
            self.s_Fn = s_Fn(Ritzel), s_Fn(Rad)
            print(f"s_Fn = {self.s_Fn}")

            # Glg D.5.06
            def d_n(idx):
                return self.geometrie.m_n * self.z_n[idx]
            self.d_n = d_n(Ritzel), d_n(Rad)
            print(f"d_n = {self.d_n}")

            # Glg D.5.07
            def d_bn(idx):
                return self.d_n[idx] * m.cos(m.radians(self.geometrie.alpha_n))
            self.d_bn = d_bn(Ritzel), d_bn(Rad)
            print(f"d_bn = {self.d_bn}")

            # Glg D.5.08
            def d_an(idx):
                return self.d_n[idx] + self.geometrie.d_a[idx] - self.geometrie.d[idx]
            self.d_an = d_an(Ritzel), d_an(Rad)
            print(f"d_an = {self.d_an}")

            # Glg D.5.09
            def alpha_an(idx):
                return m.degrees(m.acos(self.d_bn[idx] / self.d_an[idx]))
            self.alpha_an = alpha_an(Ritzel), alpha_an(Rad)
            print(f"alpha_an = {self.alpha_an}")

            # Glg D.5.10
            def y_a(idx):
                return 1 / self.z_n[idx] * (m.pi / 2 + 2 * self.geometrie.x[idx] * m.tan(m.radians(self.geometrie.alpha_n))) + involute(self.geometrie.alpha_n) - involute(self.alpha_an[idx])
            self.y_a = y_a(Ritzel), y_a(Rad)
            print(f"y_a = {self.y_a}")

            # Glg D.5.11
            def alpha_Fan(idx):
                return self.alpha_an[idx] - m.degrees(self.y_a[idx])
            self.alpha_Fan = alpha_Fan(Ritzel), alpha_Fan(Rad)
            print(f"alpha_Fan = {self.alpha_Fan}")

            # Glg D.5.12
            def h_Fa(idx):
                return self.geometrie.m_n * (0.5 * self.z_n[idx] * (m.cos(m.radians(self.geometrie.alpha_n)) / m.cos(m.radians(self.alpha_Fan[idx])) - m.cos(m.pi / 3 - m.radians(self.delta[idx])))
                                   + 0.5 * (self.geometrie.rho_fP / self.geometrie.m_n - self.G[idx] / m.cos(m.radians(self.delta[idx]))))
            self.h_Fa = h_Fa(Ritzel), h_Fa(Rad)
            print(f"h_Fa = {self.h_Fa}")

            # Glg D.5.13
            def rho_F(idx):
                return self.geometrie.rho_fP + self.geometrie.m_n * 2 * self.G[idx]**2 / m.cos(m.radians(self.delta[idx])) / (self.z_n[idx] * m.cos(m.radians(self.delta[idx]))**2 - 2 * self.G[idx])
            self.rho_F = rho_F(Ritzel), rho_F(Rad)
            print(f"rho_F = {self.rho_F}")
        else:
            raise NotImplementedError

        # Glg D.3.01
        def Y_Fa(idx):
            return (6 * self.h_Fa[idx] / self.geometrie.m_n * m.cos(m.radians(self.alpha_Fan[idx]))) / ((self.s_Fn[idx] / self.geometrie.m_n)**2 * m.cos(m.radians(self.geometrie.alpha_n)))
        self.Y_Fa = Y_Fa(Ritzel), Y_Fa(Rad)
        print(f"Y_Fa = {self.Y_Fa}")

        # Glg D.4.02
        def L_a(idx):
            return self.s_Fn[idx] / self.h_Fa[idx]
        self.L_a = L_a(Ritzel), L_a(Rad)
        print(f"L_a = {self.L_a}")

        # Glg D.4.03
        def q_s(idx):
            return self.s_Fn[idx] / 2 / self.rho_F[idx]
        self.q_s = q_s(Ritzel), q_s(Rad)
        print(f"q_s = {self.q_s}")
        assert all(1 <= q_s < 8 for q_s in self.q_s)

        # Glg D.4.01
        def Y_Sa(idx):
            return (1.2 + 0.13 * self.L_a[idx]) * m.pow(self.q_s[idx], 1 / (1.21 + 2.3 / self.L_a[idx]))
        self.Y_Sa = Y_Sa(Ritzel), Y_Sa(Rad)
        print(f"Y_Sa = {self.Y_Sa}")

        # Glg 5.08
        def Y_FS(idx):
            return self.Y_Sa[idx] * self.Y_Fa[idx]
        self.Y_FS = Y_FS(Ritzel), Y_FS(Rad)
        print(f"Y_FS = {self.Y_FS}")

        # Glg 5.10
        self.Y_beta = 1 - min(self.geometrie.epsilon_beta, 1.) * min(self.geometrie.beta, 30) / 120
        print(f"Y_beta = {self.Y_beta}")
 
        # Tabelle 5.8
        epsilon_alphan = self.geometrie.epsilon_alpha / m.cos(m.radians(self.geometrie.beta_b))**2
        def Y_S(idx):
            return self.Y_Sa[idx] * (0.6 + 0.4 * epsilon_alphan)
        self.Y_S = Y_S(Ritzel), Y_S(Rad)
        print(f"Y_S = {self.Y_S}")

        # Abschnitt 5.6
        def Y_deltarelTstat(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, ):
                raise NotImplementedError
            elif wsart in (Art.vergüteterStahl, ):
                raise NotImplementedError
            elif wsart in (Art.einsatzgehärteterStahl, ):
                return 0.44 * self.Y_S[idx] + 0.12
            elif wsart in (Art.nitrierterStahl, Art.nitrokarburierterStahl ):
                return 0.20 * self.Y_S[idx] + 0.60
            raise NotImplementedError
        # Glg 5.11, 5.12
        def Y_deltarelTdyn(idx):
            if self.q_s[idx] >= 1.5:
                return 1.
            else:
                return 0.95
        self.Y_deltarelTstat = Y_deltarelTstat(Ritzel), Y_deltarelTstat(Rad)
        self.Y_deltarelTdyn = Y_deltarelTdyn(Ritzel), Y_deltarelTdyn(Rad)
        print(f"Y_deltarelTstat = {self.Y_deltarelTstat}")
        print(f"Y_deltarelTdyn = {self.Y_deltarelTdyn}")

        # Abschnitt 5.7
        self.Y_RrelTstat = 1.
        def Y_RrelTdyn(idx):
            if self.R_z[idx] <= 16:
                return 1.
            else:
                return 0.9
        self.Y_RrelTdyn = Y_RrelTdyn(Ritzel), Y_RrelTdyn(Rad)
        print(f"Y_RrelTstat = {self.Y_RrelTstat}")
        print(f"Y_RrelTdyn = {self.Y_RrelTdyn}")

        # Tabelle 5.1
        def Y_Xdyn(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, Art.vergüteterStahl):
                if self.geometrie.m_n <= 5:
                    return 1.
                elif self.geometrie.m_n < 30:
                    return 1.03 - 0.006  * self.geometrie.m_n
                else:
                    return 0.85
            elif wsart in (Art.einsatzgehärteterStahl, Art.nitrierterStahl, Art.nitrokarburierterStahl):
                if self.geometrie.m_n <= 5:
                    return 1.
                elif self.geometrie.m_n < 25:
                    return 1.05 - 0.01  * self.geometrie.m_n
                else:
                    return 0.8
            raise NotImplementedError
        self.Y_Xstat = 1., 1.
        self.Y_Xdyn = Y_Xdyn(Ritzel), Y_Xdyn(Rad)
        print(f"Y_Xstat = {self.Y_Xstat}")
        print(f"Y_Xdyn = {self.Y_Xdyn}")

        # Tabelle 5.2
        def Y_NT(idx : int):
            wsart = self.werkstoff[idx].art
            Art = Werkstoff.Art
            if wsart in (Art.Baustahl, Art.vergüteterStahl):
                return 2.5, 1.0
            elif wsart in (Art.einsatzgehärteterStahl, ):
                return 2.5, 1.0
            elif wsart in (Art.nitrierterStahl, ):
                return 1.6, 1.0
            elif wsart in (Art.nitrokarburierterStahl, ):
                return 1.1, 1.0
            raise NotImplementedError
        Y_NTritzel = Y_NT(Ritzel)
        Y_NTrad = Y_NT(Rad)
        self.Y_NTstat = Y_NTritzel[0], Y_NTrad[0]
        self.Y_NTdyn = Y_NTritzel[1], Y_NTrad[1]

        # Glg 5.03
        def sigma_FGstat(idx):
            return self.werkstoff[idx].sigma_FE * self.Y_NTstat[idx] * self.Y_deltarelTstat[idx] * self.Y_RrelTstat * self.Y_Xstat[idx]
        def sigma_FGdyn(idx):
            return self.werkstoff[idx].sigma_FE * self.Y_NTdyn[idx] * self.Y_deltarelTdyn[idx] * self.Y_RrelTdyn[idx] * self.Y_Xdyn[idx]
        self.sigma_FGstat = sigma_FGstat(Ritzel), sigma_FGstat(Rad)
        self.sigma_FGdyn = sigma_FGdyn(Ritzel), sigma_FGdyn(Rad)
        print(f"sigma_FGstat = {self.sigma_FGstat}")
        print(f"sigma_FGdyn = {self.sigma_FGdyn}")

        # Glg 5.02
        def sigma_F0(idx):
            return self.F_t / self.geometrie.b / self.geometrie.m_n * self.Y_FS[idx] * self.Y_epsilon * self.Y_beta
        self.sigma_F0 = sigma_F0(Ritzel), sigma_F0(Rad)
        print(f"sigma_F0 = {self.sigma_F0}")

        # Glg 5.01
        def sigma_Fstat(idx):
            return self.sigma_F0[idx] * self.K_S * self.K_V[idx] * self.K_Fbeta[idx] * self.K_Falpha[idx]
        def sigma_Fdyn(idx):
            return self.sigma_F0[idx] * self.K_A * self.K_V[idx] * self.K_Fbeta[idx] * self.K_Falpha[idx]
        self.sigma_Fstat = sigma_Fstat(Ritzel), sigma_Fstat(Rad)
        self.sigma_Fdyn = sigma_Fdyn(Ritzel), sigma_Fdyn(Rad)
        print(f"sigma_Fstat = {self.sigma_Fstat}")
        print(f"sigma_Fdyn = {self.sigma_Fdyn}")

        # Glg 5.07
        def S_Fstat(idx):
            return self.sigma_FGstat[idx] / self.sigma_Fstat[idx]
        def S_Fdyn(idx):
            return self.sigma_FGdyn[idx] / self.sigma_Fdyn[idx]
        self.S_Fstat = S_Fstat(Ritzel), S_Fstat(Rad)
        self.S_Fdyn = S_Fdyn(Ritzel), S_Fdyn(Rad)
        print(f"S_Fstat = {self.S_Fstat}")
        print(f"S_Fdyn = {self.S_Fdyn}")

        print()
        return

    def ist_sicher(self, S_Hstat : float | tuple[float, float], S_Hdyn : float | tuple[float, float], S_Fstat : float | tuple[float, float], S_Fdyn : float | tuple[float, float], output = True):

        def check_value(string, val, min_or_interv):
            if isinstance(min_or_interv, tuple):
                print(f"{string} = {val} soll in {min_or_interv} liegen")
                return min_or_interv[0] <= val <= min_or_interv[1]
            else:
                print(f"{string} = {val} soll >= {min_or_interv}")
                return min_or_interv <= val

        def check_values(idx):
            result = True
            print("Ritzel" if idx == Ritzel else "Rad")
            if not check_value("S_Hstat", self.S_Hstat[idx], S_Hstat):
                result = False
                print("\033[91mstatische Grübchensicherheit ist nicht erfüllt\033[0m")
            if not check_value("S_Hdyn", self.S_Hdyn[idx], S_Hdyn):
                result = False
                print("\033[91mdynamische Grübchensicherheit ist nicht erfüllt\033[0m")
            if not check_value("S_Fstat", self.S_Fstat[idx], S_Fstat):
                result = False
                print("\033[91mstatische Zahnbruchsicherheit ist nicht erfüllt\033[0m")
            if not check_value("S_Fdyn", self.S_Fdyn[idx], S_Fdyn):
                result = False
                print("\033[91mdynamische Zahnbruchsicherheit ist nicht erfüllt\033[0m")
            print()
            return result

        with HiddenPrints(not output):
            return check_values(Ritzel) & check_values(Rad)





# getriebe = Getriebe(m_n = 3,
#                     z = (17, 38),
#                     x = (0.2, -0.2),
#                     bezugsprofil = Normalprofil2,
#                     beta = 15,
#                     k = 0,
#                     b = 40)

# getriebe.Zahnräder(P = 7.3,
#             n_1 = 1274.118,
#             verzahnungsqualität = 6,
#             werkstoff = Werkstoff(Werkstoff.Art.einsatzgehärteterStahl, 1300, 620, 220),
#             K_A = 1.1,
#             K_S = 2.1,
#             A = 0.023,

#             K_Halpha = (1, 1),
#             K_Hbeta = (1.3, 1.3),
            
#             K_Fbeta = (1, 1),
#             Z_LVRdyn = 0.85)

# result = getriebe.ist_sicher(S_Hstat_min, S_Hdyn_interval, S_Fstat_min, S_Fdyn_interval)
# sys.exit()



geometrie = Geometrie(m_n = 4,
                    z = (19, 104),
                    x = (0.5, 0.15),
                    bezugsprofil = Normalprofil1,
                    beta = 0,
                    k = 0,
                    b_d_1_verhältnis = 0.64)

assert geometrie.b / geometrie.m_n <= 30 # Konstruktionsvorgaben Tabelle 4

getriebe = Getriebe(geometrie = geometrie, P = 55,
            n_1 = 980,
            verzahnungsqualität = (6, 7),
            werkstoff = Werkstoff(Werkstoff.Art.einsatzgehärteterStahl, 1500, 860, 220),
            K_A = 1.75,
            K_S = 2.5,
            R_z = (5, 5),
            s_pr = 0,

            A = Tabelle_3_2.ohneBreitenballigkeltOderEndrücknahme,
            f_ma = (0, 0),  # Annahme, siehe Fußnote 5
            s = (0, 0),
            fertigungsverfahren = (Fertigungsverfahren.geläpptGeschliffenGeschabt, Fertigungsverfahren.geläpptGeschliffenGeschabt)
            )

result = getriebe.ist_sicher(S_Hstat_min, S_Hdyn_interval, S_Fstat_min, S_Fdyn_interval)

assert getriebe.R_z100 < 4 # Konstruktionsvorgaben Seite 7

sys.exit()




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
