# Quantum Batteries
This is an educational project. A brief description of it follows. 

# Introduction
A Quantum Battery (QB) can be defined as a $d$-dimensional quantum systems with non-degenerate energy levels from which work can be reversibly extracted – and on which energy can be reversibly deposited – by means of
cyclic unitary operations [[1]](#QB)

# Cavity Assisted Charging
One can use a [Dicke model](https://en.wikipedia.org/wiki/Dicke_model) to powerfully charge an array of 2-level systemscoupled with a quantized single-mode electro-magnetic field [[2]](#CAC). The model considered in such a case is given by the time-dependent Dicke Hamiltonian

$$
H^{(n)}=\hbar\omega_c a^\dagger a +\omega_aS_z+2\omega_c\lambda_tS_x(a+a^\dagger)
$$

where, in particular:

$$
S_i=\frac{\hbar}{2}\sum_{l=1}^n\sigma_i^{(l)}
$$
 
are the components of the collective spin operators expressed in terms of Pauli operators $\sigma_i^{(l)}$ of the $l$-th 2-level system.

# Our model
We study the following Hamiltonian:

$$
H = \omega_z S_z - g S^2_X
$$

where $S_{X, Z}=\sum_{i=1}^N\sigma^i_{X,Z}$ that can be obtained from the Dicke model [[3]](#DH1), [[4]](#DH2).

#

...

## References
<a id="QB">[1]</a> 
Francesco Campaioli, Felix A. Pollock, and Sai Vinjanampathy 2018).
*Quantum Batteries*. arXiv
[arXiv:1805.05507](https://doi.org/10.48550/arXiv.1805.05507).

<a id="CAC">[2]</a> 
Dario Ferraro, Michele Campisi, Gian Marcello Andolina, Vittorio Pellegrini, and Marco Polini (2018).
*High-Power Collective Charging of a Solid-State Quantum Battery*.
Phys. Rev. Lett. 120, 117702 
[10.1103/PhysRevLett.120.117702](https://doi.org/10.1103/PhysRevLett.120.117702).


<a id="DH1">[3]</a> 
Diego Barberena, Lucas Lamata, and Enrique Solano (2017).
*Dispersive Regimes of the Dicke Model*.
Scientific Reports volume 7, Article number: 8774 
[10.1038/s41598-017-09110-7](https://doi.org/10.1038/s41598-017-09110-7).


<a id="DH2"">[4]</a> 
Juan Román-Roche and David Zueco (2022).
*Effective theory for matter in non-perturbative cavity QED*.
SciPost Phys. Lect.Notes 50
[10.21468/SciPostPhysLectNotes.50](https://scipost.org/10.21468/SciPostPhysLectNotes.50).

