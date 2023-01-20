# Quantum Batteries
This is an educational project. A brief description follows. More info on the results in `presentation.ipynb`.

## Introduction
A Quantum Battery (QB) can be defined as a $d$-dimensional quantum systems with non-degenerate energy levels from which work can be reversibly extracted – and on which energy can be reversibly deposited – by means of cyclic unitary operations [[1]](#QB).

## Cavity Assisted Charging
One can use a [Dicke model](https://en.wikipedia.org/wiki/Dicke_model) to powerfully charge an array of 2-level systems (TLSs) coupled with a quantized single-mode electro-magnetic field [[2]](#CAC). The model considered in such a case is given by the time-dependent Dicke Hamiltonian

$$
H^{(n)}=\hbar\omega_c a^\dagger a +\omega_aS_z+2\omega_c\lambda_tS_x(a+a^\dagger)
$$

where, in particular

$$
S_i=\frac{\hbar}{2}\sum_{l=1}^n\sigma_i^{(l)}
$$
 
are the components of the collective spin operators expressed in terms of Pauli operators $\sigma_i^{(l)}$ of the $l$-th TLS.

## The model we use
We study the following Hamiltonian

$$
H = \omega_z S_z - g S^2_X
$$

that can be obtained from the Dicke model [[3]](#DH1), [[4]](#DH2).

## Final goal
We note that the state of all spin down corresponds to the fundamental state of the free Hamiltonian (without coupling) and represents an uncharged QB. We want to observe how much we can charge it (spin flip) using the evolution of the full Hamiltonian (with coupling).

## To do list
- [x] Perform the time evolution from the spin down state, implementing trotterization if needed.
- [x] Study the case $gN < 1$ and $gN >1$ (rescale for $\omega_z$). Note that the effective coupling $G=gN$ is the transition phase parameter. We should observe *universal curves* varying it.
  - [x] Study it keeping the coupling $g$ fixed.
  - [x] Study it keeping the number $N$ of TLSs fixed. 
- [x] Study the evolution of the probability of having all spin up.
- [x] Look for energy amplification: evaluate average energy of first $H$ contribution (magnetization time variation), that is $E_i(t)\equiv\text{tr}[\rho(t)H_i]-\text{tr}[\rho(0)H_i]$.
- [x] Noiseless and noisy.

### References
<a id="QB">[1]</a> 
Francesco Campaioli, Felix A. Pollock, and Sai Vinjanampathy 2018).
*Quantum Batteries*. arXiv.
[arXiv:1805.05507](https://doi.org/10.48550/arXiv.1805.05507).

<a id="CAC">[2]</a> 
Dario Ferraro, Michele Campisi, Gian Marcello Andolina, Vittorio Pellegrini, and Marco Polini (2018).
*High-Power Collective Charging of a Solid-State Quantum Battery*.
Phys. Rev. Lett. 120, 117702. 
[10.1103/PhysRevLett.120.117702](https://doi.org/10.1103/PhysRevLett.120.117702).


<a id="DH1">[3]</a> 
Diego Barberena, Lucas Lamata, and Enrique Solano (2017).
*Dispersive Regimes of the Dicke Model*.
Scientific Reports volume 7, Article number: 8774.
[10.1038/s41598-017-09110-7](https://doi.org/10.1038/s41598-017-09110-7).


<a id="DH2">[4]</a> 
Juan Román-Roche and David Zueco (2022).
*Effective theory for matter in non-perturbative cavity QED*.
SciPost Phys. Lect.Notes 50.
[10.21468/SciPostPhysLectNotes.50](https://scipost.org/10.21468/SciPostPhysLectNotes.50).


