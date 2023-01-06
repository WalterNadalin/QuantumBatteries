from numpy import abs
from qiskit.opflow import One, Zero
from itertools import combinations
from qiskit.quantum_info import SparsePauliOp
from qiskit.opflow import PauliSumOp

def U(hamiltonian: object, time: float) -> object:
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

	return exponent.exp_i()
    
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
	evolution = U(hamiltonian, time)
	evolved_state = evolution @ initial_state
	amplitude = (~state @ evolved_state).eval()

	return abs(amplitude) ** 2
	
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

def H(spins: int, frequency: float = 1, coupling: float = 1) -> object:
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
	XX = SparsePauliOp.from_sparse_list([('XX', pair, 2 * coupling) for pair in combinations(range(spins), 2)], spins)

	return PauliSumOp(Z - XX) # Return Hamiltonian
