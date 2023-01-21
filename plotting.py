from plotly.graph_objects import Figure, Contour
from tools.quantum_tools import trotter_circuit, quantum_simulator
from tools.classical_tools import classical_simulator
from tools.noise_tools import get_noise_model
from numpy import pi, linspace, arange, array, zeros_like, divide
from qiskit import transpile
from qiskit import IBMQ
from matplotlib.pyplot import title, show
from qiskit.providers.aer import QasmSimulator
from seaborn import heatmap
from matplotlib.pyplot import plot, xlabel, ylabel, title, show, legend, ylim, text, subplots, figure, grid
from matplotlib.lines import Line2D
    
def info_string(n, m, g):
    return r'$n = $' + str(n) + ', ' + r'$m = $' + str(m) + ', ' + r'$g = $' + str(g)

def compare_energies(t, qt, i_nrg, qi_nrg, c_nrg, qc_nrg, n, m, g):
    plot(t, i_nrg, '--', label = r'Classical simulation of $E_a / (n\omega_z)$')
    plot(t, c_nrg, '--', label = r'Classical simulation of $E_b / (n\omega_z)$')
    plot(qt, qi_nrg, '.', label = r'Quantum simulation of $E_a / (n\omega_z)$')
    plot(qt, qc_nrg, '.', label = r'Quantum simulation of $E_b / (n\omega_z)$')
    legend(loc = 'center right', bbox_to_anchor = (1.60, 0.5))
    title('Classical and quantum simulations of average energies \n' + info_string(n, m, g))
    xlabel(r'$\omega_z t$')
    grid()
    show()

def compare_noisy(t, qt, classical, ideal, noisy, n, m, g):
    plot(t, classical, '--', label = r'Classical simulation$')
    plot(qt, ideal, '.', label = r'Quantum ideal simulation')
    plot(qt, noisy, '.', label = r'Quantum noisy simulation')
    title('Ideal and noisy simulations of average internal energy \n' + info_string(n, m, g))
    xlabel(r'$\omega_z t$')
    grid()
    show()
    
def compare_probabilities(t, qt, p, qp, n, m, g):
    plot(t, p, '--', label = "Classical simulation")
    plot(qt, qp, '.', label = "Quantum simulation")
    xlabel(r'$\omega_z \cdot t$')
    title(r'Probability of state $|0\dots0\rangle$' + '\n' + info_string(n, m, g))
    legend(loc = 'center right', bbox_to_anchor = (1.45, 0.5))
    grid()
    show()
    
def compare_errors_4(times, quantum_times, i_nrg, n, m, g, shots, par, device_backend):
    fig, axs = subplots(2, 2, figsize = (8, 8))

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, g, noisy_backend, device_backend, shots)
    axs[0, 0].plot(quantum_times, qi_nrg, 'y.')
    axs[0, 0].plot(times, i_nrg)
    axs[0, 0].set_title('All errors\n' + info_string(n, m, g))
    axs[0, 0].grid()

    noisy_backend = QasmSimulator(noise_model = get_noise_model(par[0], 0, 0, 1, 1, n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, g, noisy_backend, device_backend, shots)
    axs[0, 1].plot(quantum_times, qi_nrg, 'g.')
    axs[0, 1].plot(times, i_nrg)
    axs[0, 1].set_title('Measurement error\n' + info_string(n, m, g))
    axs[0, 1].grid()

    noisy_backend = QasmSimulator(noise_model = get_noise_model(0, par[1], par[2], 1, 1, n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, g, noisy_backend, device_backend, shots)
    axs[1, 0].plot(quantum_times, qi_nrg, 'r.')
    axs[1, 0].plot(times, i_nrg)
    axs[1, 0].set_title('Depolarizing error')
    axs[1, 0].grid()

    noisy_backend = QasmSimulator(noise_model = get_noise_model(0, 0, 0, par[3], par[4], n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, g, noisy_backend, device_backend, shots)
    axs[1, 1].plot(quantum_times, qi_nrg, 'm.')
    axs[1, 1].plot(times, i_nrg)
    axs[1, 1].set_title('Phase Amplitude Damping')
    axs[1, 1].grid()

    axs[1, 1].set(xlabel = r'$\omega_z t$')
    axs[1, 0].set(xlabel = r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[0, 0].set(ylabel = r'$E_{a}(t)$')
    show()
    
def compare_errors_2(times, quantum_times, n, m, g, shots, par, device_backend):
    fig, axs = subplots(1, 2, figsize = (10, 10))

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[0], par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n[0], m[0], g[0], noisy_backend, device_backend, shots)
    _, i_nrg, _ = classical_simulator(times, n[0], g[0])
    axs[0].plot(quantum_times, qi_nrg, 'y.')
    axs[0].plot(times, i_nrg)
    axs[0].set_title('Noisy simulation\n' + info_string(n[0], m[0], g[0]))
    axs[0].grid()
    axs[0].set_box_aspect(1)

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[1], par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n[1], m[1], g[1], noisy_backend, device_backend, shots)
    _, i_nrg, _ = classical_simulator(times, n[1], g[1])
    axs[1].plot(quantum_times, qi_nrg, 'g.')
    axs[1].plot(times, i_nrg)
    axs[1].set_title('Noisy simulation\n' + info_string(n[1], m[1], g[1]))
    axs[1].grid()
    axs[1].set_box_aspect(1)
    
    axs[0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[1].set(xlabel=r'$\omega_z t$')
        
    show()
    
def compare_errors_our_mit(times, quantum_times, n, m, g, shots, par, device_backend):
    fig, axs = subplots(1, 2, figsize = (11, 11))
    
    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[0], par[-1]))
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[0], m[0], g[0], noisy_backend, device_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[0], m[0], g[0], noisy_backend, device_backend, shots, our_mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[0], g[0])
    axs[0].plot(quantum_times, nqi_nrg, 'y.')
    axs[0].plot(quantum_times, mnqi_nrg, 'g.')
    axs[0].plot(times, i_nrg)
    axs[0].set_title('Parity mitigation\n' + info_string(n[0], m[0], g[0]))
    axs[0].grid()
    axs[0].set_ylabel(r'$E_{a}(t)$')
    axs[0].set_box_aspect(1)

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[1], par[-1]))
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[1], m[1], g[1], noisy_backend, device_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[1], m[1], g[1], noisy_backend, device_backend, shots, our_mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[1], g[1])
    axs[1].plot(quantum_times, nqi_nrg, 'y.', label = 'Without mitigation')
    axs[1].plot(quantum_times, mnqi_nrg, 'g.', label = 'With mitigation')
    axs[1].plot(times, i_nrg)
    axs[1].set_title('Parity mitigation\n' + info_string(n[1], m[1], g[1]))
    axs[1].grid()
    axs[1].set_box_aspect(1)
    axs[1].legend(loc = 'center right', bbox_to_anchor = (1.60, 0.5))
    
    axs[0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[1].set(xlabel=r'$\omega_z t$')
        
    show()
    
def compare_errors_mit(times, quantum_times, n, m, g, shots, par, device_backend):
    fig, axs = subplots(1, 2, figsize = (11, 11))
    
    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[0], par[-1]))
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[0], m[0], g[0], noisy_backend, device_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[0], m[0], g[0], noisy_backend, device_backend, shots, measure_mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[0], g[0])
    axs[0].plot(quantum_times, nqi_nrg, 'y.')
    axs[0].plot(quantum_times, mnqi_nrg, 'g.')
    axs[0].plot(times, i_nrg)
    axs[0].set_title('Measure mitigation\n' + \
          r'$n = $' + str(n[0]) + \
          '\n' + r'$m = $' + str(m[0]) + \
          '\n' + r'$g = $' + str(g[0]))
    axs[0].grid()
    axs[0].set_ylabel(r'$E_{a}(t)$')
    axs[0].set_box_aspect(1)

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[1], par[-1]))
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[1], m[1], g[1], noisy_backend, device_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[1], m[1], g[1], noisy_backend, device_backend, shots, measure_mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[1], g[1])
    axs[1].plot(quantum_times, nqi_nrg, 'y.', label = 'Without mitigation')
    axs[1].plot(quantum_times, mnqi_nrg, 'g.', label = 'With mitigation')
    axs[1].plot(times, i_nrg)
    axs[1].set_title('Measure mitigation\n' + \
          r'$n = $' + str(n[1]) + \
          '\n' + r'$m = $' + str(m[1]) + \
          '\n' + r'$g = $' + str(g[1]))
    axs[1].grid()
    axs[1].set_box_aspect(1)
    axs[1].legend(loc = 'center right', bbox_to_anchor = (1.60, 0.5))
    
    axs[0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[1].set(xlabel=r'$\omega_z t$')
    
    show()

from qiskit.visualization import plot_histogram

def compare_histograms(n, m, g, t, parameters, shots, device_backend):
    fig, axs = subplots(1, 2, figsize = (11, 11))
    noisy_backend = QasmSimulator(noise_model = get_noise_model(*parameters[:5], n[0], parameters[-1]))
    circuit = trotter_circuit(n[0], t[0], m[0], g[0])
    circuit.measure_all()
    circuit = transpile(circuit, device_backend)
    job = noisy_backend.run(circuit, shots = shots)
    result = job.result()
    counts = result.get_counts()
    plot_histogram(counts, ax = axs[0])
    axs[0].set_box_aspect(1)
    axs[0].set_title(info_string(n[0], m[0], g[0]))
        
    circuit = trotter_circuit(n[1], t[1], m[1], g[1])
    circuit.measure_all()
    circuit = transpile(circuit, device_backend)
    job = noisy_backend.run(circuit, shots = shots)
    result = job.result()
    counts = result.get_counts()
    plot_histogram(counts, title = info_string(n[1], m[1], g[1]), ax = axs[1])
    axs[1].set_box_aspect(1)
    axs[1].set(xlabel='')
    fig.tight_layout(pad=2.0)
    show()
    
    
def compare_2_energies(t, qt, i_nrg, first, second, n, m, g, title, first_label, second_label):
    fig, axs = subplots(1, 2, figsize = (10, 10))

    axs[0].plot(t, i_nrg[0])
    axs[0].plot(qt, first[0], '.')
    axs[0].plot(qt, second[0], '.')
    axs[0].set_title(title + '\n' + info_string(n[0], m[0], g[0]))
    axs[0].grid()
    axs[0].set_box_aspect(1)

    axs[1].plot(t, i_nrg[1], label = 'Classical simulation')
    axs[1].plot(qt, first[1], '.', label = first_label)
    axs[1].plot(qt, second[1], '.', label = second_label)
    axs[1].set_title(title + '\n' + info_string(n[1], m[1], g[1]))
    axs[1].grid()
    axs[1].set_box_aspect(1)
    axs[1].legend(loc = 'center right', bbox_to_anchor = (1.60, 0.5))
    
    axs[0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[1].set(xlabel=r'$\omega_z t$')
        
    show()