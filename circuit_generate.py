import numpy as np

from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

def generate_circuit_paper(genes, qubit, parameter):
    i, q = 0, 0
    qc_r, qc_h, qc_cnot = 0, 0, 0
    qc = QuantumCircuit(qubit)
    p = ParameterVector('p', parameter)
    
    for i in range(parameter): qc.rz(0*p[i],i % qubit)   

    for gene in genes:   
        if gene[0] == 0 and gene[1] == 0 and gene[2] == 0:
            qc.h(q % qubit)
            qc_h+=1

        elif gene[0] == 0 and gene[1] == 0 and gene[2] == 1:
            qc.cx(q % qubit, (q+1) % qubit)
            qc_cnot+=1

        elif gene[0] == 0 and gene[1] == 1 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0: 
                qc.rx(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1: 
                qc.rx((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0: 
                qc.rx((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1: 
                qc.rx((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
        
        elif gene[0] == 1 and gene[1] == 1 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0: 
                qc.ry(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1: 
                qc.ry((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0: 
                qc.ry((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1: 
                qc.ry((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1

        elif gene[0] == 1 and gene[1] == 0 and gene[2] == 0:
            if gene[3] == 0 and gene[4] == 0: 
                qc.rz(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1: 
                qc.rz((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0: 
                qc.rz((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1: 
                qc.rz((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
        q+=1
    return qc, qc_r, qc_h, qc_cnot


def generate_circuit_2local(genes, qubit, parameter):
    i, q = 0, 0
    qc_r, qc_h, qc_cnot = 0, 0, 0
    qc = QuantumCircuit(qubit)
    p = ParameterVector('p', parameter)
    
    for i in range(parameter): qc.rz(0*p[i],i % qubit)   

    for gene in genes:   
        if gene[0] == 0 and gene[1] == 0 and gene[2] == 0:
            qc.h(q % qubit)
            qc_h+=1

        elif gene[0] == 0 and gene[1] == 0 and gene[2] == 1:
            qc.cx(q % qubit, (q+1) % qubit)
            qc_cnot+=1

        elif gene[0] == 1 and gene[1] == 0 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 0: 
                qc.rx(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 1: 
                qc.rx((7*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 0:
                qc.rx((3*np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 1: 
                qc.rx((5*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 0:
                qc.rx((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 1: 
                qc.rx((3*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 0: 
                qc.rx((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 1: 
                qc.rx((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
        elif gene[0] == 1 and gene[1] == 1 and gene[2] == 0:
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 0: 
                qc.ry(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 1: 
                qc.ry((7*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 0:
                qc.ry((3*np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 1: 
                qc.ry((5*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 0:
                qc.ry((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 1: 
                qc.ry((3*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 0: 
                qc.ry((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 1: 
                qc.ry((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1

        elif gene[0] == 1 and gene[1] == 1 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 0: 
                qc.rz(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 1: 
                qc.rz((7*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 0:
                qc.rz((3*np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 1: 
                qc.rz((5*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 0:
                qc.rz((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 1: 
                qc.rz((3*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 0: 
                qc.rz((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 1: 
                qc.rz((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
        q+=1
    return qc, qc_r, qc_h, qc_cnot


def generate_circuit_2local_swap(genes, qubit, parameter):
    i, q = 0, 0
    qc_r, qc_h, qc_cnot, qc_swap = 0, 0, 0, 0
    qc = QuantumCircuit(qubit)
    p = ParameterVector('p', parameter)
    
    for i in range(parameter): qc.rz(0*p[i],i % qubit)   

    for gene in genes:   
        if gene[0] == 0 and gene[1] == 0 and gene[2] == 0:
            qc.h(q % qubit)
            qc_h+=1

        elif gene[0] == 0 and gene[1] == 0 and gene[2] == 1:
            if gene[3] == 0: 
                qc.cx(q % qubit, (q+1) % qubit)
                qc_cnot+=1
            if gene[3] == 1:
                qc.swap(q % qubit, (q+1) % qubit)
                qc_swap+=1

        elif gene[0] == 1 and gene[1] == 0 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 0: 
                qc.rx(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 1: 
                qc.rx((7*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 0:
                qc.rx((3*np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 1: 
                qc.rx((5*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 0:
                qc.rx((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 1: 
                qc.rx((3*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 0: 
                qc.rx((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 1: 
                qc.rx((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
        elif gene[0] == 1 and gene[1] == 1 and gene[2] == 0:
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 0: 
                qc.ry(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 1: 
                qc.ry((7*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 0:
                qc.ry((3*np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 1: 
                qc.ry((5*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 0:
                qc.ry((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 1: 
                qc.ry((3*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 0: 
                qc.ry((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 1: 
                qc.ry((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1

        elif gene[0] == 1 and gene[1] == 1 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 0: 
                qc.rz(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 0 and gene[5] == 1: 
                qc.rz((7*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 0:
                qc.rz((3*np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 0 and gene[4] == 1 and gene[5] == 1: 
                qc.rz((5*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 0:
                qc.rz((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 0 and gene[5] == 1: 
                qc.rz((3*np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 0: 
                qc.rz((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
            if gene[3] == 1 and gene[4] == 1 and gene[5] == 1: 
                qc.rz((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_r+=1
        q+=1
    return qc, qc_r, qc_h, qc_cnot, qc_swap