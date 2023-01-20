from qiskit.providers.aer.noise import NoiseModel, pauli_error, depolarizing_error, phase_amplitude_damping_error
from numpy import exp

def get_noise_model(p_meas: float, p_dep: float, p_cnot: float, t1: float, t2: float, nqubits: int, time_gate: float):
    error_meas = pauli_error([('X',p_meas), ('I', 1 - p_meas)]) # Readout error
    depolarization = depolarizing_error(p_dep, 1)  # Single qubit depolarizing

    # Relaxation
    tph = t1 * t2 / (2 * t1 - t2)
    p1 = 1 - exp(-time_gate / t1)
    pph = 1 - exp(-time_gate / tph)
    pz = (1 - p1) * pph
    relaxation = phase_amplitude_damping_error(p1, pz, excited_state_population = 1)
    single_qubit_gate_error = depolarization.compose(relaxation) # Single qubit relaxation
    relax_cnot = relaxation.tensor(relaxation) #two qubits relaxation
    
    # Depolarization
    dep = depolarizing_error(p_cnot,1)
    dep_cnot = dep.tensor(dep) 
    two_qubits_gate_error = dep_cnot.compose(relax_cnot) # two qubits depolarizing             
    
    # Adding errors to noise model
    noise_m = NoiseModel()
    
    # Single qubit error
    for i in range(nqubits):
        noise_m.add_quantum_error(single_qubit_gate_error, ["x","sx"], [i])
        noise_m.add_quantum_error(error_meas, "measure", [i])
       
    # Two qubit error
    for j in range(nqubits):
        for k in range(nqubits):
            if (k == j + 1 or k == j - 1):         
                noise_m.add_quantum_error(two_qubits_gate_error, "cx", [j, k])
 
    return noise_m