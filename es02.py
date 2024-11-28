import numpy as np
import math
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, state_fidelity
from qiskit.providers.basic_provider import BasicProvider
from qiskit_ibm_runtime.fake_provider import FakeGuadalupeV2
from itertools import product
import pandas as pd
import os
import time
from qiskit.transpiler import CouplingMap
from qiskit.transpiler.passes import BasicSwap
from qiskit.transpiler import PassManager

# Start timer
start_time = time.time()
print(f"Starting program... {start_time}")

# Paths and files
ABS_PATH = "/home/qhd24_8/gr8_lab4/qasm_files/"
file_names = [
    "adder_small.qasm", "adder_medium.qasm", "alu-bdd_288.qasm", 
    "linearsolver.qasm", "phaseest.qasm", "rd84_253.qasm", 
    "sym10_262.qasm", "urf5_280.qasm"
]
file_names = ["adder_small.qasm"]

csv_file = '/home/qhd24_8/gr8_lab4/es01/es02-FRA.csv'


# Inizializza i risultati
results = []

# Backend simulato e gate nativi
backend = FakeGuadalupeV2()
native_gates = backend.configuration().basis_gates


#Funzione che traduce i gate non nativi in una sequenza di gate nativi
def translate_gate_to_native(gate, qargs):
    """Traduci un gate non nativo in una sequenza di gate nativi."""
    translated_gates = []
    if gate.name == "h":
            # Hadamard -> RZ(pi/2) -> SX -> RZ(pi/2)
            translated_gates.append(("rz", [math.pi / 2], qargs))
            translated_gates.append(("sx", [], qargs))
            translated_gates.append(("rz", [math.pi / 2], qargs))
    elif gate.name == "rx":
            # RX(theta) -> RZ(-pi/2) -> SX -> RZ(pi/2)
            theta = gate.params[0]
            translated_gates.append(("rz", [-math.pi / 2], qargs))
            translated_gates.append(("sx", [], qargs))
            translated_gates.append(("rz", [math.pi / 2], qargs))
            if theta != math.pi / 2:
                # Adjust for the rotation amount
                translated_gates.append(("rz", [theta], qargs))
    elif gate.name == "ry":
            # RY(theta) -> RZ(pi/2) -> SX -> RZ(theta - pi/2)
            theta = gate.params[0]
            translated_gates.append(("rz", [math.pi / 2], qargs))
            translated_gates.append(("sx", [], qargs))
            translated_gates.append(("rz", [theta - math.pi / 2], qargs))
    elif gate.name == "cz":
            # CZ -> H(target) -> CX -> H(target)
            control, target = qargs
            translated_gates.append(("rz", [math.pi / 2], [target]))
            translated_gates.append(("sx", [], [target]))
            translated_gates.append(("rz", [math.pi / 2], [target]))
            translated_gates.append(("cx", [], [control, target]))
            translated_gates.append(("rz", [math.pi / 2], [target]))
            translated_gates.append(("sx", [], [target]))
            translated_gates.append(("rz", [math.pi / 2], [target]))
    elif gate.name == "ccx":  # Toffoli (CCX)
            # Scomposizione del Toffoli (CCX) in CNOT e RZ
            control1, control2, target = qargs
            translated_gates.extend([
                ("cx", [], [control1, target]),  # CX1
                ("rz", [math.pi / 4], [target]),  # RZ(pi/4) sul target
                ("cx", [], [control2, target]),  # CX2
                ("rz", [math.pi / 4], [control2]),  # RZ(pi/4) sul control2
                ("cx", [], [control1, target]),  # CX3
                ("rz", [math.pi / 4], [target]),  # RZ(pi/4) sul target
                ("cx", [], [control2, target]),  # CX4
                ("rz", [math.pi / 4], [control2]),  # RZ(pi/4) sul control2
                ("cx", [], [control1, target])])   # CX5
    elif gate.name == "t":
            # T -> RZ(pi/4)
            translated_gates.append(("rz", [math.pi / 4], qargs))
    elif gate.name == "tdg":
            # Tdg -> RZ(-pi/4)
            translated_gates.append(("rz", [-math.pi / 4], qargs))
    elif gate.name == "u3":
            # U3(theta, phi, lambda) -> RZ(phi) -> RX(theta) -> RZ(lambda)
            theta, phi, lambda_ = gate.params
            translated_gates.append(("rz", [phi], qargs))
            translated_gates.append(("sx", [], qargs))
            translated_gates.append(("rz", [lambda_], qargs))
    else:
            raise ValueError(f"Gate {gate.name} non supportato e traduzione non definita.")
        
    return translated_gates


#Funzione che legge il file
def process_qasm_file(qasm_file, native_gates):
    """Leggi un file QASM, verifica i gate e traduci quelli non nativi."""
    # Carica il circuito dal file QASM
    circuit = QuantumCircuit.from_qasm_file(qasm_file)
    num_qubits = circuit.num_qubits      # è già un parametro di dominio pubblico, che posso chiamare anche fuori dalla funzione?
    print(num_qubits)
    # Crea un nuovo circuito con gli stessi qubit
    translated_circuit = QuantumCircuit(circuit.num_qubits)
    
    # Itera sui gate nel circuito
    for instruction in circuit.data:
        gate = instruction[0]  # Gate/istruzione
        qargs = instruction[1]  # Lista dei qubit associati
        
        if gate.name in native_gates:
            # Se il gate è nativo, aggiungilo direttamente
            translated_circuit.append(gate, qargs)
        else:
            # Se il gate non è nativo, traducilo
            translated_gates = translate_gate_to_native(gate, qargs)
            for t_gate in translated_gates:
                translated_circuit.append(t_gate, qargs)
    
    return translated_circuit

# 
def allowedConnections(backend):     # il grafo è uguale per ogni singola operazione? probabilmente i single qubit erano possibili in ogni posizione. A noi serve capire
        """Crea una semplice struttura che in futuro ci permette di controllare se due qubit sono connessi"""
#     for instr in backend.instructions:
#         print(instr)
    return connections


#Funzione che gestisce la connettività (SWAPS)
 #DA FINIRE
def routing():
        """_summary_
        """

    







# Loop sui file QASM
for name in file_names:
    filename = os.path.join(ABS_PATH, name)
    print(f"Processing file: {filename}")
    
    translated_circuit, num_qubits = process_qasm_file(filename, native_gates)
    print(translated_circuit)
    print(num_qubits)
    
    # Salvare la topologia delle connessioni. backend.instructions dizionario 
    # ---> Una lista di liste: in posizione i salvi tutti gli elementi connessi a quel qubit fisico
    connections = allowedConnections(backend)     # il grafo è uguale per ogni singola operazione? probabilmente i single qubit erano possibili in ogni posizione. A noi serve capire
    for instr in backend.instructions:
        print(instr)
    
    # Trivial mapping
    #map = [i in range()]
    
    
    
    # Prima di applicare l'istruzione "compilata", 
    # isConnected()  che ti dice se è un'operazione lecita (sono connessi i due qubit di una ): restituisce 1 (True)
    if isConneceted():   #applica istruzione sui qubit fisici 
        # cx(0,1)   ---> cx(map[0],map[1])  perché ora il qubit logico 0 è mappato nel qubit fisico map[i]

    # Altrimenti routing: SWAP(finché non avvicini) e scambiare l'array
    
    
    
    # Fidelity check: abbiamo fatto un buon lavoro?
    # Validate it using the provided qasm file, using the fidelity compared to the not compiled circuit results.

        


    


"""
    _# Memorizzare i risultati
         results.append({
        "Name of the circuit": file_name,
        "Original Depth": original_depth,
        "Original Gate Count": original_gate_count,
        "Compiled Depth": compiled_depth,
        "Compiled Gate Count": compiled_gate_count,
        "Statevector Fidelity": statevector_fidelity,
        "Probability Fidelity": probability_fidelity
    })
"""
       

# Salvare i risultati in un file CSV
df = pd.DataFrame(results)
df.to_csv(csv_file, index=False)
print(f"Results saved to {csv_file}")

# Compute total time of execution
end_time = time.time()
total_time = end_time - start_time
print(f"Total time of execution: {total_time:.2f} seconds")


    
