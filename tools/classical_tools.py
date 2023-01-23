from numpy import kron, real, abs, array
from itertools import combinations, product
from numpy import linspace, zeros_like, zeros
from scipy.linalg import expm

operator_i = array([[1, 0], [0, 1]], dtype = float) 
operator_x = array([[0, 1], [1, 0]], dtype = float)
operator_z = array([[1, 0], [0, -1]], dtype = float)

def cross_product(operators: list) -> list:
    '''
    Given a list of matrices (operators) returns the cross product between all of them as another
    matrix.
    '''
    product = kron(operators[0], operators[1]) # Cross product between the first two matrices

    for operator in operators[2:]:
        product = kron(product, operator) # Cross product between the remaining ones

    return product

def operator_single_sz(spins: int, index: int) -> list:
    '''
    Returns the spin angular momentum operator, represented as a matrix, along the z direction acting 
    on the `index` two-level system.
    '''
    single_sz = [operator_i if i != index else operator_z / 2 for i in range(spins)]
    return cross_product(single_sz)

def operator_single_sxx(spins: int, first: int, second: int) -> list:
    '''
    Returns the spin angular momentum operator, represented as a matrix, along the x direction acting 
    on the `first` and `second` two-level systems.
    '''
    single_sxx = [operator_i for _ in range(spins)] 
    single_sxx[first] = single_sxx[second] = operator_x / 2

    return cross_product(single_sxx)

def operator_sz(spins: int) -> list:
    '''
    Returns the total square spin angular momentum operator along the z direction as a matrix. 
    '''
    sz = operator_single_sz(spins, 0)
    
    for i in range(1, spins):
        sz += operator_single_sz(spins, i)
        
    return sz

def operator_sxx(spins: int) -> list:
    '''
    Returns the total square spin angular momentum operator along the x direction as a matrix. 
    '''
    sxx = 0

    for pair in combinations(range(spins), 2): # Cycling on the pairs without repetitions
        sxx += operator_single_sxx(spins, *pair)

    return sxx

def expectation_value(psi: object, operator: list) -> float:
    '''
    Returns the expected value of the observable `operator` computed on the state `psi`.
    '''
    return real(psi.conjugate() @ operator @ psi)

def state_probability(evolved_state: list, state: list) -> float:
    '''
    Computes probability of observing the state `state` given the state of the system.
    '''
    return abs(state.conjugate() @ evolved_state) ** 2

def classical_simulator(times: list, spins: int, coupling: float) -> (list, list, list):
    '''
    Given a time discretization, the number of two-level systems and the parameters of the Dicke 
    system, simulates the evolution of the system and calculate the average value of the observables
    of interest at each time step. 
    '''
    # Initial state definition
    initial_state = zeros(2 ** spins, dtype = float)
    initial_state[-1] = 1

    # State of all spins up
    state_up = zeros(2 ** spins, dtype = float)
    state_up[0] = 1

    # Dicke Hamiltonian 
    H0 = operator_sz(spins)
    H1 = -2 * coupling * operator_sxx(spins)
    hamiltonian = H0 + H1

    # Datasets containing the measures
    probabilities = zeros_like(times, dtype = float)
    internal_energy = zeros_like(times, dtype = float)
    coupling_energy = zeros_like(times, dtype = float)
    
    # Simulations
    for i, t in enumerate(times):
        evolved_state = expm(-1j * hamiltonian * t) @ initial_state
        probabilities[i] = state_probability(evolved_state, state_up)
        internal_energy[i] = expectation_value(evolved_state, H0) / spins + 1 / 2
        coupling_energy[i] = expectation_value(evolved_state, H1) / spins
        
    return probabilities, internal_energy, coupling_energy