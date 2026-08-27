import numpy as np
from scipy.special import jv, yv, jvp, yvp


# ============================================================
# FUNÇÕES DE BESSEL
# ============================================================

def Z(m, x, der=0):
    """
    Z_m(x) = J_m(x)

    der = 0 -> J_m
    der = 1 -> J'_m
    der = 2 -> J''_m
    """
    return jvp(m, x, n=der)


def W(m, x, der=0):
    """
    W_m(x) = Y_m(x)

    der = 0 -> Y_m
    der = 1 -> Y'_m
    der = 2 -> Y''_m
    """
    return yvp(m, x, n=der)


# ============================================================
# MATRIZ 8 x 8
# ============================================================

def matriz_C(omega, k, m, a, b, rho, cL, cT):
    
    mu = rho * cT**2
    lam = rho * cL**2 - 2.0 * mu

    # --------------------------------------------------------
    # Números de onda
    # --------------------------------------------------------

    kL = omega / cL
    kT = omega / cT

    # alpha_1^2 = kL^2 - k^2
    # beta_1^2  = kT^2 - k^2

    alpha1 = np.sqrt(complex(kL**2 - k**2))
    beta1  = np.sqrt(complex(kT**2 - k**2))

    # --------------------------------------------------------
    # Argumentos das funções de Bessel
    # --------------------------------------------------------

    xa = alpha1 * a
    xb = alpha1 * b

    ya = beta1 * a
    yb = beta1 * b

    C = np.zeros((8, 8), dtype=complex)

    # ========================================================
    # SUPERFÍCIE INTERNA: r = a
    # ========================================================

    # --------------------------------------------------------
    # LINHA 1: sigma_rr(a) = 0
    # --------------------------------------------------------

    C[0, 0] = (
        -lam * kL**2 * Z(m, xa)
        + 2.0 * mu * alpha1**2 * Z(m, xa, 2)
    )

    C[0, 1] = (
        -lam * kL**2 * W(m, xa)
        + 2.0 * mu * alpha1**2 * W(m, xa, 2)
    )

    C[0, 2] = (
        2.0 * mu * 1j * k * beta1 * Z(m+1, ya, 1)
    )

    C[0, 3] = (
        2.0 * mu * 1j * k * beta1 * W(m+1, ya, 1)
    )

    C[0, 4] = (
        -2.0 * mu * 1j * k * beta1 * Z(m-1, ya, 1)
    )

    C[0, 5] = (
        -2.0 * mu * 1j * k * beta1 * W(m-1, ya, 1)
    )

    C[0, 6] = (
        2.0 * mu * 1j * m * beta1 / a * Z(m, ya, 1)
        - 2.0 * mu * 1j * m / a**2 * Z(m, ya)
    )

    C[0, 7] = (
        2.0 * mu * 1j * m * beta1 / a * W(m, ya, 1)
        - 2.0 * mu * 1j * m / a**2 * W(m, ya)
    )

    # --------------------------------------------------------
    # LINHA 2: sigma_rtheta(a) = 0
    # --------------------------------------------------------

    C[1, 0] = (
        2.0 * mu * 1j * m * alpha1 / a * Z(m, xa, 1)
        - 2.0 * mu * 1j * m / a**2 * Z(m, xa)
    )

    C[1, 1] = (
        2.0 * mu * 1j * m * alpha1 / a * W(m, xa, 1)
        - 2.0 * mu * 1j * m / a**2 * W(m, xa)
    )

    C[1, 2] = (
        mu * k * beta1 * Z(m+1, ya, 1)
        - mu * k * (m+1) / a * Z(m+1, ya)
    )

    C[1, 3] = (
        mu * k * beta1 * W(m+1, ya, 1)
        - mu * k * (m+1) / a * W(m+1, ya)
    )

    C[1, 4] = (
        mu * k * beta1 * Z(m-1, ya, 1)
        + mu * k * (m-1) / a * Z(m-1, ya)
    )

    C[1, 5] = (
        mu * k * beta1 * W(m-1, ya, 1)
        + mu * k * (m-1) / a * W(m-1, ya)
    )

    C[1, 6] = (
        -mu * beta1**2 * Z(m, ya, 2)
        + mu * beta1 / a * Z(m, ya, 1)
        - mu * m**2 / a**2 * Z(m, ya)
    )

    C[1, 7] = (
        -mu * beta1**2 * W(m, ya, 2)
        + mu * beta1 / a * W(m, ya, 1)
        - mu * m**2 / a**2 * W(m, ya)
    )

    # --------------------------------------------------------
    # LINHA 3: sigma_rz(a) = 0
    # --------------------------------------------------------

    C[2, 0] = (
        2.0 * mu * 1j * k * alpha1 * Z(m, xa, 1)
    )

    C[2, 1] = (
        2.0 * mu * 1j * k * alpha1 * W(m, xa, 1)
    )

    C[2, 2] = -(
        mu * beta1**2 * Z(m+1, ya, 2)
        + mu * (m+1) / a * beta1 * Z(m+1, ya, 1)
        - mu * ((m+1) / a**2 - k**2) * Z(m+1, ya)
    )

    C[2, 3] = -(
        mu * beta1**2 * W(m+1, ya, 2)
        + mu * (m+1) / a * beta1 * W(m+1, ya, 1)
        - mu * ((m+1) / a**2 - k**2) * W(m+1, ya)
    )

    C[2, 4] = (
        mu * beta1**2 * Z(m-1, ya, 2)
        - mu * (m-1) / a * beta1 * Z(m-1, ya, 1)
        + mu * ((m-1) / a**2 + k**2) * Z(m-1, ya)
    )

    C[2, 5] = (
        mu * beta1**2 * W(m-1, ya, 2)
        - mu * (m-1) / a * beta1 * W(m-1, ya, 1)
        + mu * ((m-1) / a**2 + k**2) * W(m-1, ya)
    )

    C[2, 6] = (
        -mu * m * k / a * Z(m, ya)
    )

    C[2, 7] = (
        -mu * m * k / a * W(m, ya)
    )

    # ========================================================
    # SUPERFÍCIE EXTERNA: r = b
    # ========================================================

    # --------------------------------------------------------
    # LINHA 4: sigma_rr(b) = 0
    # --------------------------------------------------------

    C[3, 0] = (
        -lam * kL**2 * Z(m, xb)
        + 2.0 * mu * alpha1**2 * Z(m, xb, 2)
    )

    C[3, 1] = (
        -lam * kL**2 * W(m, xb)
        + 2.0 * mu * alpha1**2 * W(m, xb, 2)
    )

    C[3, 2] = (
        2.0 * mu * 1j * k * beta1 * Z(m+1, yb, 1)
    )

    C[3, 3] = (
        2.0 * mu * 1j * k * beta1 * W(m+1, yb, 1)
    )

    C[3, 4] = (
        -2.0 * mu * 1j * k * beta1 * Z(m-1, yb, 1)
    )

    C[3, 5] = (
        -2.0 * mu * 1j * k * beta1 * W(m-1, yb, 1)
    )

    C[3, 6] = (
        2.0 * mu * 1j * m * beta1 / b * Z(m, yb, 1)
        - 2.0 * mu * 1j * m / b**2 * Z(m, yb)
    )

    C[3, 7] = (
        2.0 * mu * 1j * m * beta1 / b * W(m, yb, 1)
        - 2.0 * mu * 1j * m / b**2 * W(m, yb)
    )

    # --------------------------------------------------------
    # LINHA 5: sigma_rtheta(b) = 0
    # --------------------------------------------------------

    C[4, 0] = (
        2.0 * mu * 1j * m * alpha1 / b * Z(m, xb, 1)
        - 2.0 * mu * 1j * m / b**2 * Z(m, xb)
    )

    C[4, 1] = (
        2.0 * mu * 1j * m * alpha1 / b * W(m, xb, 1)
        - 2.0 * mu * 1j * m / b**2 * W(m, xb)
    )

    C[4, 2] = (
        mu * k * beta1 * Z(m+1, yb, 1)
        - mu * k * (m+1) / b * Z(m+1, yb)
    )

    C[4, 3] = (
        mu * k * beta1 * W(m+1, yb, 1)
        - mu * k * (m+1) / b * W(m+1, yb)
    )

    C[4, 4] = (
        mu * k * beta1 * Z(m-1, yb, 1)
        + mu * k * (m-1) / b * Z(m-1, yb)
    )

    C[4, 5] = (
        mu * k * beta1 * W(m-1, yb, 1)
        + mu * k * (m-1) / b * W(m-1, yb)
    )

    C[4, 6] = (
        -mu * beta1**2 * Z(m, yb, 2)
        + mu * beta1 / b * Z(m, yb, 1)
        - mu * m**2 / b**2 * Z(m, yb)
    )

    C[4, 7] = (
        -mu * beta1**2 * W(m, yb, 2)
        + mu * beta1 / b * W(m, yb, 1)
        - mu * m**2 / b**2 * W(m, yb)
    )

    # --------------------------------------------------------
    # LINHA 6: sigma_rz(b) = 0
    # --------------------------------------------------------

    C[5, 0] = (
        2.0 * mu * 1j * k * alpha1 * Z(m, xb, 1)
    )

    C[5, 1] = (
        2.0 * mu * 1j * k * alpha1 * W(m, xb, 1)
    )

    C[5, 2] = -(
        mu * beta1**2 * Z(m+1, yb, 2)
        + mu * (m+1) / b * beta1 * Z(m+1, yb, 1)
        - mu * ((m+1) / b**2 - k**2) * Z(m+1, yb)
    )

    C[5, 3] = -(
        mu * beta1**2 * W(m+1, yb, 2)
        + mu * (m+1) / b * beta1 * W(m+1, yb, 1)
        - mu * ((m+1) / b**2 - k**2) * W(m+1, yb)
    )

    C[5, 4] = (
        mu * beta1**2 * Z(m-1, yb, 2)
        - mu * (m-1) / b * beta1 * Z(m-1, yb, 1)
        + mu * ((m-1) / b**2 + k**2) * Z(m-1, yb)
    )

    C[5, 5] = (
        mu * beta1**2 * W(m-1, yb, 2)
        - mu * (m-1) / b * beta1 * W(m-1, yb, 1)
        + mu * ((m-1) / b**2 + k**2) * W(m-1, yb)
    )

    C[5, 6] = (
        -mu * m * k / b * Z(m, yb)
    )

    C[5, 7] = (
        -mu * m * k / b * W(m, yb)
    )

    # ========================================================
    # GAUGE: r = a
    # ========================================================

    C[6, 0] = 0.0
    C[6, 1] = 0.0

    C[6, 2] = -(
        (m+1) * Z(m+1, ya)
        + a * beta1 * Z(m+1, ya, 1)
    )

    C[6, 3] = -(
        (m+1) * W(m+1, ya)
        + a * beta1 * W(m+1, ya, 1)
    )

    C[6, 4] = (
        (m-1) * Z(m-1, ya)
        - a * beta1 * Z(m-1, ya, 1)
    )

    C[6, 5] = (
        (m-1) * W(m-1, ya)
        - a * beta1 * W(m-1, ya, 1)
    )

    C[6, 6] = (
        k * a * Z(m, ya)
    )

    C[6, 7] = (
        k * a * W(m, ya)
    )

    # ========================================================
    # GAUGE: r = b
    # ========================================================

    C[7, 0] = 0.0
    C[7, 1] = 0.0

    C[7, 2] = -(
        (m+1) * Z(m+1, yb)
        + b * beta1 * Z(m+1, yb, 1)
    )

    C[7, 3] = -(
        (m+1) * W(m+1, yb)
        + b * beta1 * W(m+1, yb, 1)
    )

    C[7, 4] = (
        (m-1) * Z(m-1, yb)
        - b * beta1 * Z(m-1, yb, 1)
    )

    C[7, 5] = (
        (m-1) * W(m-1, yb)
        - b * beta1 * W(m-1, yb, 1)
    )

    C[7, 6] = (
        k * b * Z(m, yb)
    )

    C[7, 7] = (
        k * b * W(m, yb)
    )

    return C

# ============================================================
# TESTE
# ============================================================

if __name__ == "__main__":

    rho = 7800.0       # kg/m^3
    cL = 5900.0        # m/s
    cT = 3200.0        # m/s

    a = 0.01          # m
    b = 0.02          # m

    m = 0             # m=0: modo axisimétrico;

    omega = 2.0 * np.pi * 100e3   # 100 kHz
    k = 100.0                     # número de onda axial [1/m]

    C = matriz_C(
        omega=omega,
        k=k,
        m=m,
        a=a,
        b=b,
        rho=rho,
        cL=cL,
        cT=cT
    )

    print("Matriz C:")
    print(C)

    print("\nDeterminante:")
    print(np.linalg.det(C))