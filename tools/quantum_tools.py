from itertools import combinations
from qiskit import QuantumCircuit, QuantumRegister, execute, transpile
from numpy import zeros_like, array, abs, isclose, pi
from qiskit.circuit import Parameter
from qiskit import transpile
from numpy import exp
from qiskit.utils.mitigation import complete_meas_cal, CompleteMeasFitter

def quantum_simulator(times: list, spins: int, trotter_steps: int, frequency: int, coupling: int, backend: object, device_backend, shots: int, mitigation = False, our_mitigation = False, layout=[]) -> (list, list, list):
    '''
    Quantum simulation of our system given a `backend` and a number of `shots`.
    '''
    circuits = parametrized_circuits(times, spins, trotter_steps, frequency, coupling)
    circuits = transpile(circuits, device_backend, initial_layout = layout)
    job = backend.run(circuits, shots = shots)
    
    # Getting the results
    result = job.result()
    
    if mitigation:
        meas_filter = error_mitigation(backend, spins, shots).filter
        mitigated_result = meas_filter.apply(result) # Counting with error mitigation
        counts = mitigated_result.get_counts() 
    else:
        counts = result.get_counts()

    half = len(times)
    first_counts = counts[:half]
    second_counts = counts[half:]
    
    return *probability_and_internal_energy(times, first_counts, spins, shots, mitigation), \
    measure_coupling_energy(times, second_counts, coupling, spins, shots)

######################################################################################################
# NOISE MITIGATION ###################################################################################
######################################################################################################

def error_mitigation(backend: object, spins: int, shots):
    qr = QuantumRegister(spins)
    meas_calibs, state_labels = complete_meas_cal(qubit_list = range(spins), qr = qr, circlabel = 'mcal')

    # Execute the calibration circuits
    job = execute(meas_calibs, backend = backend, shots = shots)
    cal_results = job.result()

    # Calculate the calibration matrix with the noise model
    meas_fitter = CompleteMeasFitter(cal_results, state_labels, qubit_list = range(spins), circlabel = 'mcal')
    
    return meas_fitter

######################################################################################################
# USING THE OBTAINED COUNTS TO FIND THE AVERAGES ENERGIES ############################################
######################################################################################################

def measure_coupling_energy(quantum_times: list, counts: list, coupling: float, spins: int, shots: int) -> list:
    '''
    Counts the results obtained from the quantum simulation to compute the average coupling energy,
    that is the expectation value of the square total spin angular momentum along the x direction
    '''
    quantum_coupling_energy = zeros_like(quantum_times)

    # Counting to compute the total average square spin operator
    for index, count in enumerate(counts):
        for sequence in range(2 ** spins): # Iterating an all `binary sequences`
            state, state_count = state_and_count(sequence, spins, count)
            
            for first, second in combinations(range(spins), 2): # Iterating on all pairs of characters
                coefficent = -1 if state[first] == state[second] else 1
                quantum_coupling_energy[index] += coefficent * state_count

        quantum_coupling_energy[index] *= coupling / spins / shots / 2
        
    return quantum_coupling_energy

def probability_and_internal_energy(quantum_times: list, counts: list, spins: int, shots: int, mitigation = False) -> (list, list):
    '''
    Counts the results obtained from the quantum simulation to compute the average internal energy,
    that is the expectation value of the total spin angular momentum along the z direction
    '''
    quantum_internal_energy = zeros_like(quantum_times)
    quantum_probabilities = zeros_like(quantum_times)
    initial_parity = parity_operator(-spins / 2)
    
    for index, count in enumerate(counts):
        state_up = bin(0)[2:].zfill(spins) # Counting to compute tre probability of all spins up
        quantum_probabilities[index] = count.get(state_up, 0) / shots
        eliminated = 0
        counted = 0
            
        for sequence in range(2 ** spins): # Counting to compute the total average magnetization
            state, state_count = state_and_count(sequence, spins, count)
            total_spin = array([1 if char == '0' else -1 for char in list(state)]).sum() / 2

            if mitigation:
                if isclose(parity_operator(total_spin), initial_parity): # Parity must conserve
                    counted += state_count
                    quantum_internal_energy[index] += total_spin * state_count
            else:
                quantum_internal_energy[index] += total_spin * state_count
                
        factor = counted if mitigation else shots
        quantum_internal_energy[index] /= spins * factor
        quantum_internal_energy[index] += 1 / 2 
        
    return quantum_probabilities, quantum_internal_energy
    
def parity_operator(total_spin: float) -> float:
    '''
    ...
    '''
    angle = (pi * total_spin) % (2 * pi)
    
    return exp(1j * angle)
    
def state_and_count(decimal: int, spins: int, count: dict) -> (list, int):     
    '''
    ...
    '''
    state_string = bin(decimal)[2:].zfill(spins) # Converting decimal to binary sequence
    state_count = count.get(state_string, 0) # Counts related to the state defined by the sequence
    state = list(state_string) # Converting the string to a list of characters

    return state, state_count

######################################################################################################
# DEFINITION OF THE CIRCUITS #########################################################################
######################################################################################################

def get_circuit(spins: int, time_step: float, frequency: float = 1, coupling: float = 1) -> object:
    '''
    Returns a circuit implementing a single step of trotterization.
    '''
    circuit = QuantumCircuit(spins)
    circuit.rz(time_step, [i for i in range(spins)])

    for first, second in combinations(range(spins), 2):
        circuit.rxx(- coupling * time_step / frequency, first, second)
    
    return circuit

def trotter_circuit(spins: int, time: float, steps: int, frequency: float = 1, coupling: float = 1) -> object:
    '''
    Returns a circuit implementing the trotterized evolution of our system with `steps` trotter
    steps, thus composes `steps` circuits.
    '''
    time_step = time / steps
    circuit = QuantumCircuit(spins)
    circuit.x([i for i in range(spins)])
    
    for _ in range(steps - 1): # For each trotter step we compose a new single circuit
        single_circuit = get_circuit(spins, time_step, frequency, coupling)
        circuit = circuit.compose(single_circuit)
        circuit.barrier()
    
    single_circuit = get_circuit(spins, time_step, frequency, coupling)
    circuit = circuit.compose(single_circuit)

    return circuit

def parametrized_circuits(times: list, spins: int, trotter_steps: int, frequency: float, coupling: float) -> list:
    '''
    ...
    '''
    theta = Parameter('θ')
    circuit = trotter_circuit(spins, theta, trotter_steps, frequency, coupling)
    
    # Circuits for measuring probability and average internal energy
    first_circuit = circuit.copy()
    first_circuit.measure_all()
    first_circuits = [first_circuit.bind_parameters({theta: time}) for time in times]

    # Circuits for measuring coupling energy
    second_circuit = circuit.copy()
    second_circuit.h([i for i in range(spins)]) # Adding the Hadamard gates to measure along x
    second_circuit.measure_all()
    second_circuits = [second_circuit.bind_parameters({theta: time}) for time in times]

    return first_circuits + second_circuits