from qiskit import QuantumCircuit, Aer
from qiskit.opflow.expectations import PauliExpectation
from itertools import combinations
from qiskit.opflow import CircuitSampler, StateFn, I, Z

def get_circuit(spins: int, time_step: float, frequency: float = 1, coupling: float = 1) -> object:
    '''
    ...
    '''
    circuit = QuantumCircuit(spins)
    circuit.rz(2 * frequency * time_step, [i for i in range(spins)])
    
    for pair in combinations(range(spins), 2):
        first, second = pair
        circuit.rxx(-4 * coupling * time_step, first, second)
    
    return circuit

def trotter_circuit(spins: int, time: float, steps: int, frequency: float = 1, coupling: float = 1) -> object:
    '''
    ...
    '''
    time_step = time / steps
    circuit = QuantumCircuit(spins)
    
    for _ in range(steps - 1):
        single_circuit = get_circuit(spins, time_step, frequency, coupling)
        circuit = circuit.compose(single_circuit)
        circuit.barrier()
    
    single_circuit = get_circuit(spins, time_step, frequency, coupling)
    circuit = circuit.compose(single_circuit)

    return circuit

def define_operator(state: list) -> object:
    '''
    ...
    '''
    operator = (I - Z) / 2 if state[0] else (I + Z) / 2

    for component in state[1:]:
        operator = operator ^ (I - Z) / 2 if component else operator ^ (I + Z) / 2

    return StateFn(operator, is_measurement = True)