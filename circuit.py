import numpy as np

from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector

def generate_circuit(genes, qubit, parameter):
    i=0  
    q=0
    qc_one = 0
    qc_two = 0
    qc = QuantumCircuit(qubit)
    p = ParameterVector('p', parameter)
    
    for i in range(parameter): qc.rz(0*p[i],i % qubit)   

    for gene in genes:   
        if gene[0] == 0 and gene[1] == 0 and gene[2] == 0:
            qc.h(q % qubit)
            qc_one+=1

        elif gene[0] == 0 and gene[1] == 0 and gene[2] == 1:
            qc.cx(q % qubit, (q+1) % qubit)
            qc_two+=1

        elif gene[0] == 0 and gene[1] == 1 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0: 
                qc.rx(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 0 and gene[4] == 1: 
                qc.rx((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 1 and gene[4] == 0: 
                qc.rx((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 1 and gene[4] == 1: 
                qc.rx((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
        
        elif gene[0] == 1 and gene[1] == 1 and gene[2] == 1:
            if gene[3] == 0 and gene[4] == 0: 
                qc.ry(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 0 and gene[4] == 1: 
                qc.ry((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 1 and gene[4] == 0: 
                qc.ry((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 1 and gene[4] == 1: 
                qc.ry((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1

        elif gene[0] == 1 and gene[1] == 0 and gene[2] == 0:
            if gene[3] == 0 and gene[4] == 0: 
                qc.rz(np.pi*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 0 and gene[4] == 1: 
                qc.rz((np.pi/2)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 1 and gene[4] == 0: 
                qc.rz((np.pi/4)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
            if gene[3] == 1 and gene[4] == 1: 
                qc.rz((np.pi/8)*p[i % parameter], q % qubit)
                i+=1
                qc_one+=1
        q+=1
    return qc, qc_one, qc_two