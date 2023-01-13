from numpy import kron, real, abs, array
from itertools import combinations
from numpy import linspace, zeros_like, zeros
from scipy.linalg import expm

operator_i = array([[1, 0], [0, 1]], dtype = float) 
operator_x = array([[0, 1], [1, 0]], dtype = float)
operator_z = array([[1, 0], [0, -1]], dtype = float)

def cross_product(operators):
    '''
    ...
    '''
    product = kron(operators[0], operators[1])

    for operator in operators[2:]:
        product = kron(product, operator)

    return product

def operator_single_sz(spins: int, index: int) -> object:
    '''
    ...
    '''
    single_sz = [operator_i if i != index else operator_z / 2 for i in range(spins)]
    return cross_product(single_sz)

def operator_single_sxx(spins: int, first: int, second: int) -> object:
    '''
    ...
    '''
    single_sxx = [operator_i for _ in range(spins)] 
    single_sxx[first] = single_sxx[second] = operator_x / 2
    return cross_product(single_sxx)

def operator_sz(spins: int) -> object:
    '''
    ...
    '''
    sz = operator_single_sz(spins, 0)
    
    for i in range(1, spins):
        sz += operator_single_sz(spins, i)
        
    return sz

def operator_sxx(spins: int) -> object:
    '''
    ...
    '''
    sxx = 0
    
    for pair in combinations(range(spins), 2):
        sxx += operator_single_sxx(spins, *pair)

    return sxx

def expectation_value(psi: object, operator: object) -> float:
    '''
    ...
    '''
    return real(psi.conjugate() @ operator @ psi)

def state_probability(evolved_state: object, state: object) -> float:
	return abs(state.conjugate() @ evolved_state) ** 2

def classical_simulator(time_interval: list, steps: int, spins: int, frequency: float, coupling: float) -> object:
    '''
    ...
    '''
    time_interval[0] /= frequency
    time_interval[1] /= frequency
    times = linspace(*time_interval, steps)

    # Initial state definition
    initial_state = zeros(2 ** spins, dtype = float)
    initial_state[-1] = 1

    # State of all spins up
    state_up = zeros(2 ** spins, dtype = float)
    state_up[0] = 1

    # Hamiltonian operator
    H0 = frequency * operator_sz(spins)
    H1 = - 2 * coupling * operator_sxx(spins)
    hamiltonian = H0 + H1

    # Measures 
    probabilities = zeros_like(times, dtype = float)
    internal_energy = zeros_like(times, dtype = float)
    coupling_energy = zeros_like(times, dtype = float)
    
    # Simulations
    for i, t in enumerate(times):
        evolved_state = expm(-1j * hamiltonian * t) @ initial_state
        probabilities[i] = state_probability(evolved_state, state_up)
        internal_energy[i] = expectation_value(evolved_state, H0) / frequency / spins + 1 / 2
        coupling_energy[i] = expectation_value(evolved_state, H1) / frequency / spins
        
    return times, probabilities, internal_energy, coupling_energy