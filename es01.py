import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector
from qiskit.providers.basic_provider import BasicProvider
from qiskit_ibm_runtime.fake_provider import FakeGuadalupeV2
from itertools import product
import pandas as pd
import os
import time

# Inizio del timer
start_time = time.time()
print(f"Starting program... {start_time}")

# Percorsi e file
ABS_PATH = "/home/qhd24_8/gr8_lab4/qasm_files/"
file_names = [
    "adder_small.qasm", "adder_medium.qasm", "alu-bdd_288.qasm", 
    "linearsolver.qasm", "phaseest.qasm", "rd84_253.qasm", 
    "sym10_262.qasm", "urf5_280.qasm"
]


# Intestazione per il file CSV
header = ["Name of the circuit", "Non compiled Circuit Depth", "Non compiled Gate Count",
          "Compiled Circuit Depth", "Compiled Gate Count", 
          "Statevector Fidelity", "Probability Fidelity"]

# Inizializza i risultati
circuit_results = []

# Loop sui file QASM
for name in file_names:
    filename = os.path.join(ABS_PATH, name)
    print(f"Processing file: {filename}")
    
    # 1) import the quantum circuit
    qc = QuantumCircuit.from_qasm_file(filename)
    
    # 2) simulate it with Statevector
    statevector = Statevector(qc)
    print(np.abs(statevector.inner(statevector)))

    # 3) copy the circuit, add the measurement and simulate with the BasicProvider;
    qc_m = qc.copy()         # copy the circuit
    qc_m.measure_all()       # add the measurement

    backend = BasicProvider().get_backend("basic_simulator")
    result = backend.run(qc_m).result()
    counts = result.get_counts()
    print("Counts:", counts)
    
    num_qubits = qc.num_qubits
    # Generate all possible binary outcomes for num_qubits
    all_states = [''.join(state) for state in product('01', repeat=num_qubits)]

    # Calculate probabilities, including zero counts
    total_shots = sum(counts.values())
    probabilities = [counts.get(state, 0) / total_shots for state in all_states]

    # 4) evaluate the circuit depth and gate count
    print("Not compiled solutions:")
    nonCompDepth = qc_m.depth()
    nonCompGateCount = qc_m.count_ops()
    print(f"Circuit depth: {nonCompDepth}")
    print(f"Gate count: {nonCompGateCount}")

    # 5) compile the circuit considering as basis gates: ry, rx, rz, cx, considering as optimization levels [0, 1, 2, 3] and compute:
    optimization_levels = [0, 1, 2, 3]
    for level in optimization_levels:
        qc_t = transpile(qc, basis_gates=['ry', 'rx', 'rz', 'cx'], optimization_level=level)
    
    print("Compiled results")
    depth = qc_t.depth()
    gateCount = qc_t.count_ops()
    print(f"Circuit depth: {depth}")
    print(f"Gate count: {gateCount}")

    # 6) compile the circuit considering the FakeGuadalupeV2, considering as optimization levels [0, 1, 2, 3] and compute
    print("Backend FakeGuadalupeV2")
    backend2 = FakeGuadalupeV2()

    # To see backend properties
    coupling_map = backend2.configuration().coupling_map
    print("Coupling map:", coupling_map)

    properties = backend2.properties()
    print("Properties:", properties)

    # Extract calibration data for each qubit
    print("\n--- Qubit Properties ---")
    for qubit_index in range(len(properties.qubits)):
        print(f"Qubit {qubit_index}:")
        for item in properties.qubits[qubit_index]:
            print(f" {item.name}: {item.value} {item.unit}")

    # Extract gate errors for each gate
    print("\n--- Gate Properties ---")
    for gate in properties.gates:
        print(f"Gate: {gate.gate} on qubits {gate.qubits}")
        for param in gate.parameters:
            print(f" {param.name}: {param.value} {param.unit}")

    #native_gates = backend2.configuration().basis_gates
    #print("Native gates (basis gates):", native_gates)

    statevector2 = Statevector(qc_t)

    qc_t.measure_all()
    result2 = backend2.run(qc_t).result()
    counts2 = result2.get_counts()

    num_qubits2 = qc_t.num_qubits
    # Generate all possible binary outcomes for num_qubits
    all_states2 = [''.join(state) for state in product('01', repeat=num_qubits2)]
    # Calculate probabilities, including zero counts
    total_shots2 = sum(counts2.values())
    print(total_shots2)
    probabilities2 = [counts2.get(state, 0) / total_shots2 for state in all_states2]
    
    # To compute fidelity between two statevectors
    fidelity = np.abs(statevector.inner(statevector2))
    print(f"Fidelity between non compiled and compiled statevectors: {fidelity} \n")
    prob_fidelity = (np.sum(np.sqrt(probabilities) * np.sqrt(probabilities2)))**2
    print(f"Fidelity between non compiled and compiled probabilities: {prob_fidelity} \n")

    # Save results in the list
    circuit_results.append([name, level, nonCompDepth, nonCompGateCount, depth, gateCount, fidelity, prob_fidelity])
# Converte i risultati in un DataFrame Pandas
df = pd.DataFrame(circuit_results, columns=header)

# Scrivi i risultati nel file CSV in modalità append
csv_file = '/home/qhd24_8/gr8_lab4/es01Gio/es01.csv'
df.to_csv(csv_file, mode='a', header=not os.path.exists(csv_file), index=False)

print("Dati salvati con successo nel file es01.csv")


# Calcola il tempo totale di esecuzione
end_time = time.time()
total_time = end_time - start_time
print(f"Tempo totale di esecuzione: {total_time:.2f} secondi")
