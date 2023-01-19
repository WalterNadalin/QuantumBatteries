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

def plotcontour(x, y, u, size = (700, 700)):
    a, b = size
    fig = Figure(data = Contour(z = u, x = x, y = y))
    fig.update_layout({
    "title": {"text": r"$\text{Maximum power } \frac{P_{max}}{N\cdot\omega_z^2}$", "x": 0.5, "y": 0.9, "font": {"size": 14} },
    "showlegend": True,
    "xaxis": {"title": "$g$", "showticklabels": True, "dtick": 0.2}, 
    "yaxis": {"title": "$N$", "showticklabels": True, "dtick": 1}, 
    "autosize": False, 
    "width":a, 
    "height":b})
    fig.show()
    
def compare_energies(times, quantum_times, internal_energy, quantum_internal_energy, coupling_energy, quantum_coupling_energy, n, spins, g):
    plot(times, internal_energy, '--', label = r'Classical simulation of $E_a / (n\omega_z)$')
    plot(times, coupling_energy, '--', label = r'Classical simulation of $E_b / (n\omega_z)$')
    plot(quantum_times, quantum_internal_energy, '.', label = r'Quantum simulation of $E_a / (n\omega_z)$')
    plot(quantum_times, quantum_coupling_energy, '.', label = r'Quantum simulation of $E_b / (n\omega_z)$')
    #styles = [':', '--']
    #lines = [Line2D([0], [0], color = 'black', linestyle = s) for s in styles]
    #labels = ['Quantum simulation', 'Classical simulation']
    legend(loc = 'center right', bbox_to_anchor = (1.60, 0.5))
    title('Classical and quantum simulations of average energies \n' + \
          r'$n = $' + str(spins) + \
          '\n' + r'$m = $' + str(n) + \
          '\n' + r'$g = $' + str(g))
    xlabel(r'$\omega_z t$')
    grid()
    show()

def compare_noisy(times, quantum_times, classical, first, second, n, spins, g):
    plot(times, classical, '--', label = r'Classical simulation$')
    plot(quantum_times, first, '.', label = r'Quantum ideal simulation')
    plot(quantum_times, second, '.', label = r'Quantum noisy simulation')
    #styles = [':', '--']
    #lines = [Line2D([0], [0], color = 'black', linestyle = s) for s in styles]
    #labels = ['Quantum simulation', 'Classical simulation']
    legend(loc = 'center right', bbox_to_anchor = (1.65, 0.5))
    title('Ideal and noisy simulations of average internal energy \n' + \
          r'$n = $' + str(spins) + \
          '\n' + r'$m = $' + str(n) + \
          '\n' + r'$g = $' + str(g))
    xlabel(r'$\omega_z t$')
    grid()
    show()

    
def compare_probabilities(times, quantum_times, probabilities, quantum_probabilities, n, spins, g):
    plot(times, probabilities, '--', label = "Classical simulation")
    plot(quantum_times, quantum_probabilities, '.', label = "Quantum simulation")
    xlabel(r'$\omega_z \cdot t$')
    title(r'Probability of state $|0\dots0\rangle$' + '\n' + \
          r'$n = $' + str(spins) + \
          '\n' + r'$m = $' + str(n) + \
          '\n' + r'$g = $' + str(g))
    legend(loc = 'center right', bbox_to_anchor = (1.45, 0.5))
    grid()
    show()
    
def compare_errors4(times, quantum_times, i_nrg, n, m, frequency, g, shots, par):
    fig, axs = subplots(2, 2, figsize = (7, 7))

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, frequency, g, noisy_backend, shots)
    axs[0, 0].plot(quantum_times, qi_nrg, 'y.')
    axs[0,0].plot(times, i_nrg)
    axs[0, 0].set_title('All errors')
    axs[0, 0].grid()

    noisy_backend = QasmSimulator(noise_model = get_noise_model(par[0], 0, 0, 1, 1, n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, frequency, g, noisy_backend, shots)
    axs[0, 1].plot(quantum_times, qi_nrg, 'g.')
    axs[0,1].plot(times, i_nrg)
    axs[0, 1].set_title('Measurement error')
    axs[0, 1].grid()

    noisy_backend = QasmSimulator(noise_model = get_noise_model(0, par[1], par[2], 1, 1, n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, frequency, g, noisy_backend, shots)
    axs[1, 0].plot(quantum_times, qi_nrg, 'r.')
    axs[1,0].plot(times, i_nrg)
    axs[1, 0].set_title('Depolarizing error')
    axs[1, 0].grid()

    noisy_backend = QasmSimulator(noise_model = get_noise_model(0, 0, 0, par[3], par[4], n, par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n, m, frequency, g, noisy_backend, shots)
    axs[1, 1].plot(quantum_times, qi_nrg, 'm.')
    axs[1,1].plot(times, i_nrg)
    axs[1, 1].set_title('Phase Amplitude Damping')
    axs[1, 1].grid()

    axs[1, 1].set(xlabel=r'$\omega_z t$')
    axs[1, 0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[0, 0].set(ylabel=r'$E_{a}(t)$')

    show()
    
def compare_errors2(times, quantum_times, n, m, frequency, g, shots, par):
    fig, axs = subplots(1, 2, figsize = (10, 10))

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[0], par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n[0], m[0], frequency, g[0], noisy_backend, shots)
    _, i_nrg, _ = classical_simulator(times, n[0], frequency, g[0])
    axs[0].plot(quantum_times, qi_nrg, 'y.')
    axs[0].plot(times, i_nrg)
    axs[0].set_title('Noisy simulation\n' + \
          r'$n = $' + str(n[0]) + \
          '\n' + r'$m = $' + str(m[0]) + \
          '\n' + r'$g = $' + str(g[0]))
    axs[0].grid()
    axs[0].set_box_aspect(1)

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[1], par[-1]))
    _, qi_nrg, _ = quantum_simulator(quantum_times, n[1], m[1], frequency, g[1], noisy_backend, shots)
    _, i_nrg, _ = classical_simulator(times, n[1], frequency, g[1])
    axs[1].plot(quantum_times, qi_nrg, 'g.')
    axs[1].plot(times, i_nrg)
    axs[1].set_title('Noisy simulation\n' + \
          r'$n = $' + str(n[1]) + \
          '\n' + r'$m = $' + str(m[1]) + \
          '\n' + r'$g = $' + str(g[1]))
    axs[1].grid()
    axs[1].set_box_aspect(1)
    
    axs[0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[1].set(xlabel=r'$\omega_z t$')
        
    show()
    
def compare_errors_our_mit(times, quantum_times, n, m, frequency, g, shots, par):
    fig, axs = subplots(1, 2, figsize = (11, 11))
    
    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[0], par[-1]))
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[0], m[0], frequency, g[0], noisy_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[0], m[0], frequency, g[0], noisy_backend, shots, our_mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[0], frequency, g[0])
    axs[0].plot(quantum_times, nqi_nrg, 'y.')
    axs[0].plot(quantum_times, mnqi_nrg, 'g.')
    axs[0].plot(times, i_nrg)
    axs[0].set_title('Parity mitigation\n' + \
          r'$n = $' + str(n[0]) + \
          '\n' + r'$m = $' + str(m[0]) + \
          '\n' + r'$g = $' + str(g[0]))
    axs[0].grid()
    axs[0].set_ylabel(r'$E_{a}(t)$')
    axs[0].set_box_aspect(1)

    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[1], par[-1]))
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[1], m[1], frequency, g[1], noisy_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[1], m[1], frequency, g[1], noisy_backend, shots, our_mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[1], frequency, g[1])
    axs[1].plot(quantum_times, nqi_nrg, 'y.', label = 'Without mitigation')
    axs[1].plot(quantum_times, mnqi_nrg, 'g.', label = 'With mitigation')
    axs[1].plot(times, i_nrg)
    axs[1].set_title('Parity mitigation\n' + \
          r'$n = $' + str(n[1]) + \
          '\n' + r'$m = $' + str(m[1]) + \
          '\n' + r'$g = $' + str(g[1]))
    axs[1].grid()
    #axs[1].set_ylabel(r'$E_{a}(t)$')
    axs[1].set_box_aspect(1)
    axs[1].legend(loc = 'center right', bbox_to_anchor = (1.60, 0.5))
    
    axs[0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[1].set(xlabel=r'$\omega_z t$')
        
    
    show()
    
def compare_errors_mit(times, quantum_times, n, m, frequency, g, shots, par):
    fig, axs = subplots(1, 2, figsize = (11, 11))
    
    noisy_backend = QasmSimulator(noise_model = get_noise_model(*par[:5], n[0], par[-1]))
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[0], m[0], frequency, g[0], noisy_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[0], m[0], frequency, g[0], noisy_backend, shots, mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[0], frequency, g[0])
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
    _, nqi_nrg, _ = quantum_simulator(quantum_times, n[1], m[1], frequency, g[1], noisy_backend, shots)
    _, mnqi_nrg, _ = quantum_simulator(quantum_times,n[1], m[1], frequency, g[1], noisy_backend, shots, mitigation = True)
    _, i_nrg, _ = classical_simulator(times, n[1], frequency, g[1])
    axs[1].plot(quantum_times, nqi_nrg, 'y.', label = 'Without mitigation')
    axs[1].plot(quantum_times, mnqi_nrg, 'g.', label = 'With mitigation')
    axs[1].plot(times, i_nrg)
    axs[1].set_title('Measure mitigation\n' + \
          r'$n = $' + str(n[1]) + \
          '\n' + r'$m = $' + str(m[1]) + \
          '\n' + r'$g = $' + str(g[1]))
    axs[1].grid()
    #axs[1].set_ylabel(r'$E_{a}(t)$')
    axs[1].set_box_aspect(1)
    axs[1].legend(loc = 'center right', bbox_to_anchor = (1.60, 0.5))
    
    axs[0].set(xlabel=r'$\omega_z t$', ylabel=r'$E_{a}(t)$')
    axs[1].set(xlabel=r'$\omega_z t$')
    
    show()
