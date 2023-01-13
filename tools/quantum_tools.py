from qiskit.opflow import I, X, Z
from itertools import combinations
from qiskit import QuantumCircuit
from qiskit.opflow.expectations import PauliExpectation
from qiskit.opflow import CircuitSampler, StateFn
from qiskit.providers.aer import QasmSimulator

expectation = PauliExpectation()

def quantum_operator_single_sz(spins: int, index: int) -> object:
    '''
    ...
    '''
    single_sz = [I if i != index else Z / 2 for i in range(spins)]
    product = quantum_cross_product(single_sz)
    return StateFn(product, is_measurement = True)

def quantum_operator_single_sxx(spins: int, first: int, second: int) -> object:
    '''
    ...
    '''
    single_sxx = [I for _ in range(spins)] 
    single_sxx[first] = single_sxx[second] = X / 2
    product = quantum_cross_product(single_sxx)
    return StateFn(product, is_measurement = True)

def measure_operator(operator: object, circuit: object, shots: int) -> object:
    '''
    ...
    '''
    backend = QasmSimulator(shots = shots)
    sampler = CircuitSampler(backend) 
    measure = expectation.convert(operator @ StateFn(circuit))
    sample = sampler.convert(measure) 
    return sample.eval().real

def state_operator(state: list) -> object:
    '''
    ...
    '''
    operator = [(I - Z) / 2 if single_state else (I + Z) / 2 for single_state in state]
    product = quantum_cross_product(operator)

    return StateFn(product, is_measurement = True)

def get_circuit(spins: int, time_step: float, frequency: float = 1, coupling: float = 1) -> object:
    '''
    ...
    '''
    circuit = QuantumCircuit(spins)
    circuit.rz(frequency * time_step, [i for i in range(spins)])
    
    for pair in combinations(range(spins), 2):
        first, second = pair
        circuit.rxx(- coupling * time_step, first, second)
    
    return circuit

def trotter_circuit(spins: int, time: float, steps: int, frequency: float = 1, coupling: float = 1) -> object:
    '''
    ...
    '''
    time_step = time / steps
    circuit = QuantumCircuit(spins)
    circuit.x([i for i in range(spins)])
    
    for _ in range(steps - 1):
        single_circuit = get_circuit(spins, time_step, frequency, coupling)
        circuit = circuit.compose(single_circuit)
        circuit.barrier()
    
    single_circuit = get_circuit(spins, time_step, frequency, coupling)
    circuit = circuit.compose(single_circuit)

    return circuit

def quantum_cross_product(operators):
    product = operators[0] ^ operators[1]

    for operator in operators[2:]:
        product = product ^ operator

    return product