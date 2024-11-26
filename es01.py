import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector
from qiskit.providers.basic_provider import BasicProvider
from qiskit_ibm_runtime.fake_provider import FakeGuadalupeV2
from itertools import product
import pandas as pd
import os
import time

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


# Header for the CSV file
header = ["Name of the circuit", "Non compiled Circuit Depth", "Non compiled Gate Count",
          "Compiled Circuit Depth", "Compiled Gate Count", 
          "Statevector Fidelity", "Probability Fidelity"]

# Initialize results
circuit_results = []

# Define the two backends
backend = BasicProvider().get_backend("basic_simulator")
backend2 = FakeGuadalupeV2()

# Loop through QASM files
for name in file_names:
    filename = os.path.join(ABS_PATH, name)
    print(f"Processing file: {filename}")
    
    # 1) import the quantum circuit
    qc = QuantumCircuit.from_qasm_file(filename)
    
    # 2) simulate it with Statevector
    statevector = Statevector(qc)
    print(f"Initial statevector: {statevector}")

    # 3) copy the circuit, add the measurement and simulate with the BasicProvider;
    qc_m = qc.copy()         # copy the circuit
    qc_m.measure_all()       # add the measurement

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
    print(f"Circuit depth (non-compiled): {nonCompDepth}")
    print(f"Gate count (non-compiled): {nonCompGateCount}")

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

    # To see backend properties
    coupling_map = backend2.configuration().coupling_map
    # print("Coupling map:", coupling_map)

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

# Converts results into a Pandas DataFrame
df = pd.DataFrame(circuit_results, columns=header)

# Write the results to the CSV file in "append" mode
csv_file = '/home/qhd24_8/gr8_lab4/es01Gio/es01.csv'
df.to_csv(csv_file, mode='a', header=not os.path.exists(csv_file), index=False)

print("Data successfully saved in the file es01.csv")


# Compute total time of execution
end_time = time.time()
total_time = end_time - start_time
print(f"Total time of execution: {total_time:.2f} seconds")
