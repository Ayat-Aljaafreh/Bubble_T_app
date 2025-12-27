import math


def main():

    n = int(input("inter the number of components: ")) #number of component


    while True:
        nonideal_liq = input("Is liquid phase non-ideal? (y/n): ").strip().lower()
        if nonideal_liq in ["y", "n"]:
            nonideal_liq = (nonideal_liq == "y")
            break
        print("Enter y or n.")

    while True:
        nonideal_gas = input("Is gas phase non-ideal? (y/n): ").strip().lower()
        if nonideal_gas in ["y", "n"]:
            nonideal_gas = (nonideal_gas == "y")
            break
        print("Enter y or n.")

    if nonideal_liq:
        a = [[0.0 for j in range(n)] for i in range(n)]
        print("\nEnter Wilson equation constants a_ij:")
        for i in range(n):
            for j in range(n):
                if i != j:
                    a[i][j] = float(input(f"a{i+1}{j+1}: "))

        V = []                                             #molar volume
        print("\nEnter molar volumes V_i:")
        for i in range(n):
            V.append(float(input(f"V{i+1}: ")))
        R = float(input("Enter gas constant R: "))
    else:
        a = [[0.0 for j in range(n)] for i in range(n)]
        V = [1.0] * n
        R = float(1)



    if nonideal_gas:
        R = float(input("Enter gas constant R: "))

        Bik = [[0.0 for j in range(n)] for i in range(n)]
        print("\nEnter second virial coefficients Bik (including Bii):")
        for i in range(n):
            for k in range(i, n):  # الجزء العلوي فقط بما فيها القطر
                value = float(input(f"B{i+1}{k+1}: "))
                Bik[i][k] = value
                if i != k:
                    Bik[k][i] = value  # لضمان التماثل
    else:
        Bik = [[0.0 for j in range(n)] for i in range(n)]









    A = []                                             #Antoine constants
    B = []
    C = []

    print("\nEnter Antoine constants:")
    for i in range(n):
        A.append(float(input(f"A{i+1}: ")))
        B.append(float(input(f"B{i+1}: ")))
        C.append(float(input(f"C{i+1}: ")))

    Ptot = float(input("\nEnter system pressure P: "))


    ALLOWED_T_UNITS = ["c", "k"]                      #Antoine Equation Units
    ALLOWED_P_UNITS = ["mmhg", "bar", "kpa"]
    while True:
        T_unit_Ant = input("Enter temperature unit for Antoine equation (C or K): ").strip().lower()

        if T_unit_Ant in ALLOWED_T_UNITS:
            break
        else:
            print("Invalid unit. Allowed temperature units: C, K.")
    while True:
        P_unit_Ant = input("Enter saturation pressure unit (mmHg, bar, kPa): ").strip().lower()

        if P_unit_Ant in ALLOWED_P_UNITS:
            break
        else:
            print("Invalid unit. Allowed pressure units: mmHg, bar, kPa.")

    P = convert_pressure_to_kpa(Ptot, P_unit_Ant)

    while True:


        # Ask user which component to fix
        while True:
            try:
                k_index = int(input(f"\nEnter the number of the fixed component (1-{n}): "))
                if 1 <= k_index <= n:
                    k_index -= 1
                    break
                else:
                    print(f"Invalid input. Enter a number between 1 and {n}.")
            except ValueError:
                print("Invalid input. Enter an integer.")



        while True:                                      #mole fractions
            x = []
            print("\nEnter mole fractions x_i:")

            for i in range(n):
                x.append(float(input(f"x{i+1}: ")))

            if abs(sum(x))== 1:
                break
            else:
                print(f"Error: Sum of mole fractions = {sum(x):.6f}")
                print("Please re-enter mole fractions so that their sum equals 1.")



            # Ask user for error tolerance in percentage
        while True:
            try:
                tol_percent = float(input("\nEnter error tolerance (%) for iteration (e.g., 0.01 for 0.01%): "))
                if tol_percent > 0:
                    break
                else:
                    print("Tolerance must be positive.")
            except ValueError:
                print("Invalid input. Enter a numeric value.")

        tol = tol_percent
        print("\nSystem non-ideality options:")



        T0 =calculate_T0 (P, A, B, C, P_unit_Ant, T_unit_Ant,x)
        print(f"\nT0 = {T0:.8f} {T_unit_Ant.upper()}")





        Results = bubble_point_iteration(T0, P, x, A, B, C, T_unit_Ant, V, a, R, k_index, P_unit_Ant, tol, nonideal_liq, nonideal_gas, Bik)

        choice = input("\nEnter (1) to run again with new P, x, tolerance, or (0) to exit: ")

        if choice == "0":
            print("Program terminated.")
            break





def calculate_T0(P, A, B, C, P_unit_Ant, T_unit_Ant,x):
    Ti_sat = []                                     #Calculate Tsat
    n = len(A)

    for i in range(n):
        T_i = B[i] / (A[i] - math.log(P)) - C[i]  # in same unit as Antoine
        Ti_sat.append(T_i)

    T0 = 0.0
    for i in range(n):
        T0 += x[i] * Ti_sat[i]

    return T0


def calculate_Psat(T_1, T_unit_Ant, A, B, C, P_unit_Ant):
    Psat = []
    n = len(A)

    for i in range(n):

        Psat_i = math.exp((A[i] - B[i] / (T_1 + C[i])))
        Psat.append(Psat_i)

    return Psat



def calculate_Lambda(a, V, R, T_Ant, T_unit):

    T_K = to_kelvin(T_Ant, T_unit)
    n = len(V)
    Lambda = [[0.0 for j in range(n)] for i in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                Lambda[i][j] = 1.0
            else:
                Lambda[i][j] = (V[j] / V[i]) * math.exp(-a[i][j] / (R * T_K))

    return Lambda



def calculate_gamma(x, Lambda):
    n = len(x)
    gamma = []

    for i in range(n):

        # first sum: sum_j (xj * Lambda_ij)
        sum1 = 0.0
        for j in range(n):
            sum1 += x[j] * Lambda[i][j]

        # second term: sum_k [ xk * Lambda_ki / sum_j(xj * Lambda_kj) ]
        sum2 = 0.0
        for k in range(n):
            denom = 0.0
            for j in range(n):
                denom += x[j] * Lambda[k][j]

            sum2 += x[k] * Lambda[k][i] / denom

        ln_gamma_i = 1.0 - math.log(sum1) - sum2
        gamma.append(math.exp(ln_gamma_i))

    return gamma


def convert_pressure_to_kpa(P, unit):
    if unit == "kpa":
        return P
    elif unit == "bar":
        return P * 100.0
    elif unit == "mmhg":
        return P * 0.133322
    else:
        raise ValueError("Unsupported pressure unit")

def to_kelvin(T, unit):
    if unit == "k":
        return T
    elif unit == "c":
        return T + 273
    else:
        raise ValueError("Invalid temperature unit")

def calculate_phi(T_old, P, Bik, R, yi, T_unit):

    T_K = to_kelvin(T_old, T_unit)

    n = len(yi)
    phi = []

    # Compute Sij = 2*Bik - Bii - Bkk
    S = [[2 * Bik[i][j] - Bik[i][i] - Bik[j][j] for j in range(n)] for i in range(n)]

    for k in range(n):
        sum_term = 0.0
        for i in range(n):
            for j in range(n):
                sum_term += yi[i] * yi[j] * (2 * S[i][k] - S[i][j])

        ln_phi_k = (P / (R * T_K)) * (Bik[k][k] + 0.5 * sum_term)
        phi_k = math.exp(ln_phi_k)
        phi.append(phi_k)

    return phi





def calculate_alpha(Pi_sat, k_index):
    alpha = [Pi_sat[i] / Pi_sat[k_index] for i in range(len(Pi_sat))]
    return alpha



def bubble_point_iteration(T0, P, x, A, B, C, T_unit_Ant, V, a, R, k_index, P_unit_Ant, tol, nonideal_liq, nonideal_gas, Bik, max_iter=100):



        # Pi_sat for each component using latest temperature (Tk)


    T_old = T0  # <-- initial guess = T0
    iteration = 0
    phi = [1.0] * len(x)
    Pi_sat = calculate_Psat(T_old, T_unit_Ant, A, B, C, P_unit_Ant)
    if nonideal_liq:
        Lambda = calculate_Lambda(a, V, R, T_old, T_unit_Ant)
        gamma = calculate_gamma(x, Lambda)
    else:
        gamma = [1.0] * len(x)

    alpha = calculate_alpha(Pi_sat, k_index)

    numerator_sum = sum(((x[i] * gamma[i] * alpha[i])/ phi[i]) for i in range(len(x)))
    Pk_sat = P / numerator_sum

    Tk_sat = (B[k_index] / (A[k_index] - math.log(Pk_sat))) - C[k_index]


    while iteration < max_iter:

        iteration += 1
        T_old = Tk_sat

        Pi_sat = calculate_Psat(T_old, T_unit_Ant, A, B, C, P_unit_Ant)

        # ---- Liquid phase ----
        # Lambda and gamma
        if nonideal_liq:
            Lambda = calculate_Lambda(a, V, R, T_old, T_unit_Ant)
            gamma = calculate_gamma(x, Lambda)
        else:
            gamma = [1.0] * len(x)

        yi = [(x[i] * gamma[i] * Pi_sat[i])/(P * phi[i]) for i in range(len(x))]

        # ---- Gas phase ----
        if nonideal_gas:
            phi = calculate_phi(T_old, P, Bik, R, yi, T_unit_Ant)
        else:
            phi = [1.0] * len(x)



        alpha = calculate_alpha(Pi_sat, k_index)

        # Pk_sat from formula
        numerator_sum = sum(((x[i] * gamma[i] * alpha[i])/ phi[i]) for i in range(len(x)))
        Pk_sat = P / numerator_sum

        # Tk_sat for the fixed component
        Tk_sat = B[k_index] / (A[k_index] - math.log(Pk_sat)) - C[k_index]


        # Calculate error (%)
        error = (abs(T_old - Tk_sat) / T_old) * 100



        print(f"\nIteration {iteration}")
        print(f"T_old = {T_old:.8f} {T_unit_Ant.upper()}")
        print(f"Tk_sat = {Tk_sat:.8f} {T_unit_Ant.upper()}")
        print(f"Error  = {error:.8f} %")
        print(f"Pk_sat = {Pk_sat:.8f} kPa")

        print("\nComponent-wise properties:")
        print("i   Pi_sat (kPa)     gamma_i        alpha_ik")
        print("-" * 60)

        for i in range(len(x)):
            print(f"{i+1:<3} {Pi_sat[i]:<14.8f} {gamma[i]:<14.8f} {alpha[i]:<14.8f}")






        # Check convergence
        if error <= tol:
            T_final = Tk_sat
            break

        # Update T_old for next iteration
        T_old = Tk_sat


    else:
        print("Warning: Maximum iterations reached without convergence.")
        T_final = Tk_sat

    # Compute yi mole fractions in vapor

    Pi_sat = calculate_Psat(T_final, T_unit_Ant, A, B, C, P_unit_Ant)

    if nonideal_liq:
        Lambda = calculate_Lambda(a, V, R, T_final, T_unit_Ant)
        gamma = calculate_gamma(x, Lambda)
    else:
        gamma = [1.0] * len(x)

    yi = [(x[i] * gamma[i] * Pi_sat[i])/(P * phi[i]) for i in range(len(x))]

    print("\nConverged Bubble Point Results:")
    print(f"Final T = {T_final:.8f} {T_unit_Ant.upper()}, Iterations = {iteration}")
    print("Vapor phase mole fractions yi:")
    for i, y in enumerate(yi):
        print(f"y{i+1} = {y:.8f}")

    return T_final, yi



if __name__ == "__main__":
    main()





