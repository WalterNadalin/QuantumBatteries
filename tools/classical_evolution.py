from numpy import kron, real, abs, array
from itertools import combinations

operator_i = array([[1, 0], [0, 1]], dtype = complex) 
operator_x = array([[0, 1], [1, 0]], dtype = complex)
operator_z = array([[1, 0], [0, -1]], dtype = complex)

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
	'''
	...
	''' 
	amplitude = state.conjugate() @ evolved_state

	return abs(amplitude) ** 2