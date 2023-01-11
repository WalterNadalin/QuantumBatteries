from numpy import abs, kron, array
from qiskit.opflow import One, Zero
from itertools import combinations
from qiskit.quantum_info import SparsePauliOp
from qiskit.opflow import PauliSumOp

operator_i = array([[1, 0], [0, 1]], dtype = complex) 
operator_x = array([[0, 1], [1, 0]], dtype = complex) 
operator_z = array([[1, 0], [0, -1]], dtype = complex) 

def evolution(hamiltonian: object, time: float) -> object:
	'''
	Returns the unitary evolution operator given the Hamiltonian.

	Parameters
	----------
	hamiltonian : Hamiltonian matrix representation.
	time : evolution time.
	  
	Returns
	-------
	Unitary evolution operator.
	'''
	exponent = hamiltonian * time

	return exponent.exp_i() # Return the exponential of `-i` times `exponent`
    
def state_probability(hamiltonian: object, time: float, initial_state: object, state: object) -> float:
	'''
	Returns the probability of observing a given state at a given time.

	Parameters
	----------
	hamiltonian : Hamiltonian matrix representation.
	time : time at which to compute the probability.
	initial_state : system initial state.
	state: state of which compute the probability.

	Returns
	-------
	Probability.
	''' 
	evolved_state = evolution(hamiltonian, time) @ initial_state
	amplitude = (~state @ evolved_state).eval() # `_state` gives the complex conjugate

	return abs(amplitude) ** 2 # `eval()` returns the inner product 
	
def define_state(state: list) -> object:
	'''
	Given a list of zeros and ones, return the state vector based on it.

	Parameters
	----------
	state : list of zeros and ones.

	Returns
	-------
	State vector.
	'''
	vector = One if state[0] else Zero

	for component in state[1:]:
		vector = vector ^ One if component else vector ^ Zero

	return vector

def dicke_hamiltonian(spins: int, frequency: float = 1, coupling: float = 1) -> object:
	'''
	Returns the Dicke Hamiltonian matrix representation that we are studying.

	Parameters
	----------
	spins : number of 2-level system.
	frequency : frequency associated to the energy gap of the 2-level system.
	coupling : intensity of the coupling with the cavity electro-magnetic field.

	Returns
	-------
	Dicke Hamiltonian matrix representation.
	'''   
	Z = SparsePauliOp.from_sparse_list([('Z', [i], frequency) for i in range(spins)], spins)
	XX = SparsePauliOp.from_sparse_list([('XX', pair, 2 * coupling) for pair in \
                                         combinations(range(spins), 2)], spins)

	return PauliSumOp(Z - XX) # Return Hamiltonian

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
    single_sz = [operator_i if i != index else operator_z for i in range(spins)]
    return cross_product(single_sz)

def operator_sz(spins: int) -> object:
    '''
    ...
    '''
    sz = operator_single_sz(spins, 0)
    
    for i in range(1, spins):
        sz += operator_single_sz(spins, i)
        
    return sz

def operator_single_sxx(spins: int, first: int, second: int) -> object:
    '''
    ...
    '''
    single_sxx = [operator_i for _ in range(spins)] 
    single_sxx[first] = single_sxx[second] = operator_x
    return cross_product(single_sxx)

def operator_sxx(spins: int) -> object:
    '''
    ...
    '''
    sxx = 0
    
    for pair in combinations(range(spins), 2):
        sxx += 2 * operator_single_sxx(spins, *pair)

    return sxx