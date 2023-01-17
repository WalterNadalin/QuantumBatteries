from itertools import combinations
from qiskit import QuantumCircuit, QuantumRegister, execute, transpile
from numpy import zeros_like, array, abs
from qiskit.circuit import Parameter
from qiskit import transpile
from qiskit.providers.aer.noise import NoiseModel, pauli_error, depolarizing_error, phase_amplitude_damping_error
from numpy import exp
from qiskit.utils.mitigation import complete_meas_cal, CompleteMeasFitter

def get_noise_model(p_meas: float, p_dep: float, p_cnot: float, t1: float, t2: float, nqubits: int, time_gate: float):
    # Readout error
    error_meas = pauli_error([('X',p_meas), ('I', 1 - p_meas)])
    
    # Single qubit depolarizing
    depolarization = depolarizing_error(p_dep, 1)

    # Single qubit relaxation
    tph = t1 * t2 / (2 * t1 - t2)
    p1 = 1 - exp(-time_gate / t1)
    pph = 1 - exp(-time_gate / tph)
    pz = (1 - p1) * pph
    relaxation = phase_amplitude_damping_error(p1, pz)
    single_qubit_gate_error = depolarization.compose(relaxation)
                    
    #two qubits relaxation
    relax_cnot = relaxation.tensor(relaxation)
             
    #two qubits depolarizing               
    dep = depolarizing_error(p_cnot,1)
    dep_cnot = dep.tensor(dep) 
    two_qubits_gate_error = dep_cnot.compose(relax_cnot)         
    
    # Adding errors to noise model
    noise_m = NoiseModel()
    
    #single qubit error
    for i in range(nqubits):
        noise_m.add_quantum_error(single_qubit_gate_error, ["x","sx"], [i])
        noise_m.add_quantum_error(error_meas, "measure", [i])
    
    # two qubit error
    for j in range(nqubits):
        for k in range(nqubits):
            if (k == j+1 or k == j-1):         
                noise_m.add_quantum_error(two_qubits_gate_error, "cx", [j, k])
    
    return noise_m

def error_mitigation(backend: object, spins: int, shots):
    qr = QuantumRegister(spins)
    meas_calibs, state_labels = complete_meas_cal(qubit_list = range(spins), qr = qr, circlabel='mcal')

    # Execute the calibration circuits
    job = execute(meas_calibs, backend = backend, shots = shots)
    cal_results = job.result()

    # Calculate the calibration matrix with the noise model
    meas_fitter = CompleteMeasFitter(cal_results, state_labels, qubit_list=range(spins), circlabel='mcal')
    
    return meas_fitter

def quantum_simulator(times: list, spins: int, trotter_steps: int, frequency: int, coupling: int, backend: object, shots: int, mitigation = False) -> (list, list, list):
    '''
    Quantum simulation of our system given a `backend` and a number of `shots`.
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

    # Transpiling and sending the job
    circuits = first_circuits + second_circuits
    circuits = transpile(circuits, backend)
    job = backend.run(circuits, shots = shots)
    
    # Getting the results
    result = job.result()
    counts = result.get_counts()
    half = len(times)
    first_counts = counts[:half]
    second_counts = counts[half:]
    
    # Extracting the wanted information from the data
    quantum_probabilities, quantum_internal_energy = \
    probability_and_internal_energy(times, first_counts, spins, shots)
    quantum_coupling_energy = measure_coupling_energy(times, second_counts, coupling, spins, shots)
    
    if mitigation:
        # Apply error mitigation
        meas_filter = error_mitigation(backend, spins, shots).filter
        mitigated_result = meas_filter.apply(result)
        mitigated_counts = mitigated_result.get_counts()
        first_mitigated_counts = mitigated_counts[:half]
        second_mitigated_counts = mitigated_counts[half:]
        
        # Counting with error mitigation
        quantum_probabilities, quantum_internal_energy = \
        probability_and_internal_energy(times, first_mitigated_counts, spins, shots)
        quantum_coupling_energy = \
        measure_coupling_energy(times, second_mitigated_counts, coupling, spins, shots)
    
    return quantum_probabilities, quantum_internal_energy, quantum_coupling_energy

def measure_coupling_energy(quantum_times: list, counts: list, coupling: float, spins: int, shots: int) -> list:
    '''
    Counts the results obtained from the quantum simulation to compute the average coupling energy,
    that is the expectation value of the square total spin angular momentum along the x direction
    '''
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
    Counts the results obtained from the quantum simulation to compute the average internal energy,
    that is the expectation value of the total spin angular momentum along the z direction
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
    Returns a circuit implementing a single step of trotterization.
    '''
    circuit = QuantumCircuit(spins)
    circuit.rz(time_step, [i for i in range(spins)])
    
    for pair in combinations(range(spins), 2):
        first, second = pair
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
    
    for _ in range(steps - 1):
        single_circuit = get_circuit(spins, time_step, frequency, coupling)
        circuit = circuit.compose(single_circuit)
        circuit.barrier()
    
    single_circuit = get_circuit(spins, time_step, frequency, coupling)
    circuit = circuit.compose(single_circuit)

    return circuit