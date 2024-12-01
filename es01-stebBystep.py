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

# Define the header for the DataFrame
header = ["Name of the circuit", "Backend", "Optimization Level", "Depth", "Gate Count", "Statevector Fidelity", "Probability Fidelity"]

# Create a DataFrame to store the results
df = pd.DataFrame(columns=header)

# Define the file path for the CSV
csv_file = '/home/qhd24_8/gr8_lab4/es01Gio/es01-stepBystep.csv'

# Verifica se il file esiste ed è vuoto
if not os.path.exists(csv_file) or os.path.getsize(csv_file) == 0:
    # Se il file non esiste o è vuoto, scrivi l'header
    df.to_csv(csv_file, mode='w', header=True, index=False)
else:
    # Se il file esiste ed ha già dei dati, scrivi solo i dati senza l'header
    df.to_csv(csv_file, mode='a', header=False, index=False)


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
    nonCompDepth = qc.depth()
    nonCompGateCount = sum(qc.count_ops().values())
    print(f"Circuit depth (non-compiled): {nonCompDepth}")
    print(f"Gate count (non-compiled): {nonCompGateCount}")
    
    # Add the non-compiled results to the DataFrame
    #df.loc[len(df)] = [name, "Basic Simulator", 'non compiled', nonCompDepth, nonCompGateCount, None, None]
    row = [name, "Basic Simulator", 'non compiled', nonCompDepth, nonCompGateCount, None, None]
    pd.DataFrame([row], columns=header).to_csv(csv_file, mode='a', header=False, index=False)

    # 5) compile the circuit considering as basis gates: ry, rx, rz, cx, considering as optimization levels [0, 1, 2, 3] and compute:
    optimization_levels = [0, 1, 2, 3]
    for level in optimization_levels:
        qc_t = transpile(qc, basis_gates=['ry', 'rx', 'rz', 'cx'], optimization_level=level)
    
        print(f"Compiled results. Level: {level}")
        depth = qc_t.depth()
        gateCount = sum(qc_t.count_ops().values())
        print(f"Circuit depth: {depth}")
        print(f"Gate count: {gateCount}")
        
        # Simulate it with Statevector (ideal simulation)
        statevector2 = Statevector(qc_t)

        # To compute fidelity between two statevectors
        statevector_fidelity = np.abs(statevector.inner(statevector2))
        print(f"Fidelity between non compiled and compiled statevectors: {statevector_fidelity} \n")
        
        # Simulate with BasicProvider (Measurement-based simulation)
        qc_tm = qc_t.copy()         # copy the circuit
        qc_tm.measure_all()       # add the measurement
    
        result2 = backend2.run(qc_tm).result()
        counts2 = result2.get_counts()

        num_qubits2 = qc_t.num_qubits
        # Generate all possible binary outcomes for num_qubits
        all_states2 = [''.join(state) for state in product('01', repeat=num_qubits2)]
        # Calculate probabilities, including zero counts
        total_shots2 = sum(counts2.values())
        print(total_shots2)
        probabilities2 = [counts2.get(state, 0) / total_shots2 for state in all_states]
        
        prob_fidelity = (np.sum(np.sqrt(probabilities) * np.sqrt(probabilities2)))**2
        print(f"Fidelity between non compiled and compiled probabilities: {prob_fidelity} \n")

        # Add results to the DataFrame
        #df.loc[len(df)] = [name, 'Basic Simulator', level, depth, gateCount, statevector_fidelity, prob_fidelity] 
        row = [name, 'Basic Simulator', level, depth, gateCount, statevector_fidelity, prob_fidelity] 
        pd.DataFrame([row], columns=header).to_csv(csv_file, mode='a', header=False, index=False)     
        

    # 6) compile the circuit considering the FakeGuadalupeV2, considering as optimization levels [0, 1, 2, 3] and compute
    print("Backend FakeGuadalupeV2")

    for level in optimization_levels:
        qc_t = transpile(qc, backend=backend, optimization_level=level)

        print(f"Compiled results. Level: {level}")
        depth = qc_t.depth()
        gateCount = sum(qc_t.count_ops().values())
        print(f"Circuit depth: {depth}")
        print(f"Gate count: {gateCount}")
        
        # Simulate it with Statevector
        statevector2 = Statevector(qc_t)
                
        # print("Stavector: \n", statevector)
        # print("Stavector compiled: \n", statevector2)
        # Simulate with BasicProvider (Measurement-based simulation)
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

        # Add results to the DataFrame
        #df.loc[len(df)] = [name, 'FakeGuadalupe', level, depth, gateCount, fidelity, prob_fidelity]
        row = [name, 'FakeGuadalupe', level, depth, gateCount, fidelity, prob_fidelity]
        pd.DataFrame([row], columns=header).to_csv(csv_file, mode='a', header=False, index=False)   
     


print("Data successfully saved in the file es01.csv")
# Compute total time of execution
end_time = time.time()
total_time = end_time - start_time
print(f"Total time of execution: {total_time:.2f} seconds")
