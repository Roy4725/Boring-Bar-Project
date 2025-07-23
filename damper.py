import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, fsolve
import math

# Parameters (you may need to adjust a and b based on your specific problem)
a = 37  # Inner radius or dimension
b = 40  # Outer radius or dimension
mu = 0.004
target_K = 115

def calculate_S(L, a, b):
    """Calculate S parameter"""
    return (2/3) * (L / (b - a)) * ((b**2 + a*b + a**2) / (b + a)**2)

def calculate_beta_approx(S, a, b):
    """Calculate beta_approx parameter"""
    numerator = 4 * math.pi
    denominator = (((S**2 + 2.45) / (S**2 + 1.45)) * math.log(b/a) - 
                   ((b**2 - a**2) / (b**2 + a**2)))
    return numerator / denominator

def calculate_K(L, a, b, mu):
    """Calculate K for given L"""
    S = calculate_S(L, a, b)
    beta_approx = calculate_beta_approx(S, a, b)
    K = beta_approx * mu * L
    return K

def objective_function(L):
    """Objective function to minimize - difference between calculated K and target K"""
    calculated_K = calculate_K(L, a, b, mu)
    return abs(calculated_K - target_K)

def K_equation(L):
    """Equation to solve: K(L) - target_K = 0"""
    return calculate_K(L, a, b, mu) - target_K

# Method 1: Using optimization to minimize the difference
print("Method 1: Optimization approach")
print(f"Target K: {target_K}")
print(f"Parameters: a = {a}, b = {b}, mu = {mu}")
print("-" * 50)

# Find optimal L using optimization
result = minimize_scalar(objective_function, bounds=(0.1, 100), method='bounded')
optimal_L = result.x
optimal_K = calculate_K(optimal_L, a, b, mu)

print(f"Optimal L: {optimal_L:.6f}")
print(f"Achieved K: {optimal_K:.6f}")
print(f"Error: {abs(optimal_K - target_K):.6f}")

# Method 2: Using root finding
print("\nMethod 2: Root finding approach")
try:
    # Try different initial guesses
    for initial_guess in [1, 5, 10, 20]:
        try:
            L_solution = fsolve(K_equation, initial_guess)[0]
            if L_solution > 0:  # Valid solution
                K_check = calculate_K(L_solution, a, b, mu)
                print(f"L solution: {L_solution:.6f}")
                print(f"K check: {K_check:.6f}")
                print(f"Error: {abs(K_check - target_K):.6f}")
                break
        except:
            continue
except Exception as e:
    print(f"Root finding failed: {e}")

# Create a graph showing K vs L relationship
L_range = np.linspace(0.1, 50, 1000)
K_values = []

for L in L_range:
    try:
        K = calculate_K(L, a, b, mu)
        K_values.append(K)
    except:
        K_values.append(np.nan)

plt.figure(figsize=(10, 6))
plt.plot(L_range, K_values, 'b-', linewidth=2, label='K(L)')
plt.axhline(y=target_K, color='r', linestyle='--', linewidth=2, label=f'Target K = {target_K}')
plt.axvline(x=optimal_L, color='g', linestyle=':', linewidth=2, label=f'Optimal L = {optimal_L:.3f}')
plt.xlabel('L')
plt.ylabel('K')
plt.title('K vs L Relationship')
plt.grid(True, alpha=0.3)
plt.legend()
plt.xlim(0, 50)
plt.ylim(0, 200)

# Add annotation for the solution point
plt.annotate(f'Solution\nL = {optimal_L:.3f}\nK = {optimal_K:.1f}', 
             xy=(optimal_L, optimal_K), 
             xytext=(optimal_L + 5, optimal_K + 20),
             arrowprops=dict(arrowstyle='->', color='red'),
             fontsize=10,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))

plt.tight_layout()
plt.show()

# Additional analysis
print(f"\nAdditional Analysis:")
print(f"S parameter at optimal L: {calculate_S(optimal_L, a, b):.6f}")
print(f"Beta_approx at optimal L: {calculate_beta_approx(calculate_S(optimal_L, a, b), a, b):.6f}")

# Test different L values around the optimal solution
print(f"\nSensitivity analysis around optimal L:")
test_L_values = [optimal_L * 0.9, optimal_L * 0.95, optimal_L, optimal_L * 1.05, optimal_L * 1.1]
for test_L in test_L_values:
    test_K = calculate_K(test_L, a, b, mu)
    print(f"L = {test_L:.4f}, K = {test_K:.2f}, Error = {abs(test_K - target_K):.2f}")