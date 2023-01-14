from itertools import combinations
from qiskit import QuantumCircuit
from numpy import zeros_like, array, abs
from qiskit.circuit import Parameter
from qiskit import transpile

def quantum_simulator(times: list, spins: int, trotter_steps: int, frequency: int, coupling: int, backend: object, shots: int) -> (list, list, list):
    theta = Parameter('θ')
    circuit = trotter_circuit(spins, theta, trotter_steps, frequency, coupling)
    
    # Circuits for measuring probability and average internal energy
    first_circuit = circuit.copy()
    first_circuit.measure_all()
    first_circuits = [first_circuit.bind_parameters({theta: time}) for time in times]

    # Circuits for measuring coupling energy
    second_circuit = circuit.copy()
    second_circuit.h([i for i in range(spins)])
    second_circuit.measure_all()
    second_circuits = [second_circuit.bind_parameters({theta: time}) for time in times]
    
    # transpiling
    first_circuits = transpile(first_circuits, backend)
    second_circuits = transpile(second_circuits, backend)
    
    # Simulating
    first_job = backend.run(first_circuits, shots = shots)
    first_counts = first_job.result().get_counts()
    second_job = backend.run(second_circuits, shots = shots)
    second_counts = second_job.result().get_counts()
    
    # Counting
    quantum_probabilities, quantum_internal_energy = \
    probability_and_internal_energy(times, first_counts, spins, shots)
    quantum_coupling_energy = \
    measure_coupling_energy(times, second_counts, coupling, spins, shots)
    
    return quantum_probabilities, quantum_internal_energy, quantum_coupling_energy

def measure_coupling_energy(quantum_times: list, counts: list, coupling: float, spins: int, shots: int) -> list:
    quantum_coupling_energy = zeros_like(quantum_times)

    # Counting to compute the total average square spin operator
    for index, count in enumerate(counts):
        for sequence in range(2 ** spins): # Iterating an all binary sequences
            state_string = bin(sequence)[2:].zfill(spins) # Converting binary sequence to a string
            state_count = count.get(state_string, 0)
            state = list(state_string) # Converting the string to a list of characters

            for pair in combinations(range(spins), 2): # Iterating on all pairs of characters
                first, second = pair
                coefficent = -1 if state[first] == state[second] else 1
                quantum_coupling_energy[index] += coefficent * state_count

        quantum_coupling_energy[index] *= coupling / spins / shots / 2
        
    return quantum_coupling_energy

def probability_and_internal_energy(quantum_times: list, counts: list, spins: int, shots: int) -> (list, list):
    '''
    ...
    '''
    quantum_internal_energy = zeros_like(quantum_times)
    quantum_probabilities = zeros_like(quantum_times)

    for index, count in enumerate(counts):
        # Counting to compute tre probability of the state: all spins up
        state_up = bin(0)[2:].zfill(spins)
        quantum_probabilities[index] = count.get(state_up, 0) / shots

        # Counting to compute the total average magnetization: same as above
        for sequence in range(2 ** spins):
            state_string = bin(sequence)[2:].zfill(spins)
            state_count = count.get(state_string, 0)
            state = list(state_string)
            coefficents = array([1 if char == '0' else -1 for char in list(state)])

            for k in range(spins):
                quantum_internal_energy[index] += coefficents[k] * state_count / shots / 2

        quantum_internal_energy[index] /= spins
        quantum_internal_energy[index] += 1 / 2 

    return quantum_probabilities, quantum_internal_energy

def get_circuit(spins: int, time_step: float, frequency: float = 1, coupling: float = 1) -> object:
    '''
    ...
    '''
    circuit = QuantumCircuit(spins)
    circuit.rz(time_step, [i for i in range(spins)])
    
    for pair in combinations(range(spins), 2):
        first, second = pair
        circuit.rxx(- coupling * time_step / frequency, first, second)
    
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