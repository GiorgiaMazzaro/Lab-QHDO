import numpy as np
import math
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt
from qiskit.quantum_info import Statevector
from qiskit.providers.basic_provider import BasicProvider
from qiskit_ibm_runtime.fake_provider import FakeGuadalupeV2
from itertools import product
import pandas as pd
import os
import time
from qiskit.transpiler import CouplingMap
from qiskit.transpiler.passes import BasicSwap
from qiskit.transpiler import PassManager
from collections import deque
from qiskit.visualization import plot_histogram

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
#file_names = ["adder_small.qasm"]

csv_file = '/home/qhd24_8/gr8_lab4/es01Gio/es02-new.csv'

# Creazione del DataFrame iniziale
header = [
    "Name of the circuit",
    "Original Depth",
    "Original Gate Count",
    "Compiled Depth",
    "Compiled Gate Count",
    "Additional swaps (x routing)",
    "Probability Fidelity",
    "Statevector Fidelity",
]
df = pd.DataFrame(columns=header)

# Backend simulato e gate nativi
backend_basic = BasicProvider().get_backend("basic_simulator")
backend = FakeGuadalupeV2()
native_gates = backend.configuration().basis_gates
max_qubits = backend.configuration().num_qubits    # Ottieni il numero di qubit disponibili dal backend


# Funzione che traduce i gate non nativi in una sequenza di gate nativi
def translate_gate_to_native(gate, qargs, qc):     # gli passiamo il circuito così li applica direttamente 
    """Traduci un gate non nativo in una sequenza di gate nativi."""
    # Determina il numero di qubit (assumiamo che tutti i qubit siano nel circuito)
    num_qubits = len(qargs)
 
    if gate.name == "h":
        # Hadamard -> RZ(pi/2) -> SX -> RZ(pi/2)
        qc.rz(math.pi / 2, qargs[0])  # Applica RZ(pi/2) su qargs[0]
        qc.sx(qargs[0])  # Applica SX su qargs[0]
        qc.rz(math.pi / 2, qargs[0])  # Applica RZ(pi/2) su qargs[0]
    elif gate.name == "rx":
        # RX(theta) -> RZ(-pi/2) -> SX -> RZ(pi/2)
        theta = gate.params[0]
        qc.rz(-math.pi / 2, qargs[0])  # Applica RZ(-pi/2) su qargs[0]
        qc.sx(qargs[0])  # Applica SX su qargs[0]
        qc.rz(math.pi / 2, qargs[0])  # Applica RZ(pi/2) su qargs[0]
        if theta != math.pi / 2:
            # Adjust for the rotation amount
            qc.rz(theta, qargs[0])  # Applica il valore di theta su qargs[0]
    elif gate.name == "ry":
        # RY(theta) -> RZ(pi/2) -> SX -> RZ(theta - pi/2)
        theta = gate.params[0]
        qc.rz(math.pi / 2, qargs[0])  # Applica RZ(pi/2) su qargs[0]
        qc.sx(qargs[0])  # Applica SX su qargs[0]
        qc.rz(theta - math.pi / 2, qargs[0])  # Applica RZ(theta - pi/2) su qargs[0]
    elif gate.name == "cz":
        # CZ -> H(target) -> CX -> H(target)
        control, target = qargs
        qc.h(target)  # Applica H su target
        qc.cx(control, target)  # Applica CX su (control, target)
        qc.h(target)  # Applica H su target
        #translated_gates.append(qc)  # Aggiungi il circuito all'elenco
    elif gate.name == "ccx":  # Toffoli (CCX)
        # Scomposizione del Toffoli (CCX) in CNOT e RZ
        control1, control2, target = qargs
        qc.cx(control1, target)  # CX1
        qc.rz(math.pi / 4, target)  # RZ(pi/4) sul target
        qc.cx(control2, target)  # CX2
        qc.rz(math.pi / 4, control2)  # RZ(pi/4) sul control2
        qc.cx(control1, target)  # CX3
        qc.rz(math.pi / 4, target)  # RZ(pi/4) sul target
        qc.cx(control2, target)  # CX4
        qc.rz(math.pi / 4, control2)  # RZ(pi/4) sul control2
        qc.cx(control1, target)  # CX5
    elif gate.name == "s":    ## AGGIUNTO
        # S -> RZ(pi/2)   ???
        qc.rz(math.pi / 2, qargs[0])  # Applica RZ(pi/2) su qargs[0]   ## DA CONTROLLARE
    elif gate.name == "t":
        # T -> RZ(pi/4)
        qc.rz(math.pi / 4, qargs[0])  # Applica RZ(pi/4) su qargs[0]
    elif gate.name == "tdg":
        # Tdg -> RZ(-pi/4)
        qc.rz(-math.pi / 4, qargs[0])  # Applica RZ(-pi/4) su qargs[0]
    elif gate.name == "u3":
        # U3(theta, phi, lambda) -> RZ(phi) -> SX -> RZ(lambda)
        theta, phi, lambda_ = gate.params
        qc.rz(phi, qargs[0])  # Applica RZ(phi) su qargs[0]
        qc.sx(qargs[0])  # Applica SX su qargs[0]
        qc.rz(lambda_, qargs[0])  # Applica RZ(lambda_) su qargs[0]
    else:
        raise ValueError(f"Gate {gate.name} non supportato e traduzione non definita.")
    
    return 0

#Funzione che legge il file
def process_qasm_file(qasm_file, native_gates):   # il numero massimo di qubit è fissato dal backend (16 per fakeguadalupe)
    """Leggi un file QASM, verifica i gate e traduci quelli non nativi."""
    # Carica il circuito dal file QASM
    circuit = QuantumCircuit.from_qasm_file(qasm_file)
    num_qubits = circuit.num_qubits      # è già un parametro di dominio pubblico, che posso chiamare anche fuori dalla funzione?
    
    # Crea un nuovo circuito con il numero di qubit del backend
    translated_circuit = QuantumCircuit(max_qubits)

    # Itera sui gate nel circuito
    for instruction in circuit.data:
        gate = instruction.operation  # Access the gate operation
        qargs = instruction.qubits  # Access the qubits    # è del tipo (Qubit(QuantumRegister(4, 'q'), 0), Qubit(QuantumRegister(4, 'q'), 1))
        # Estrai gli indici dei qubit da qargs
        qargs_indices = [q._index for q in qargs]
        
        if gate.name in native_gates:
            # Se il gate è nativo, aggiungilo direttamente
            translated_circuit.append(gate, qargs_indices)
        else:
            # Se il gate non è nativo, traducilo
            check = translate_gate_to_native(gate, qargs_indices, translated_circuit)  # check = 0 se è andato tutto bene

    return circuit, translated_circuit, num_qubits

def allowedConnections(backend):
    """Crea una lista di liste per rappresentare le connessioni tra i qubit fisici."""
    coupling_map = backend.configuration().coupling_map  # Ottieni la mappa delle connessioni

    # Inizializza la lista delle connessioni con set per rimuovere i duplicati automaticamente
    connections = [set() for _ in range(max_qubits)]

    for edge in coupling_map:
        q1, q2 = edge
        connections[q1].add(q2)
        connections[q2].add(q1)  # Connessione bidirezionale

    # Converti i set in liste per la compatibilità con il resto del codice
    return [list(conn) for conn in connections]

def isConnected(connections, q1, q2):
    """Controlla se due qubit fisici sono direttamente connessi."""
    return q2 in connections[q1]    # 1 se sono connessi

def apply_gate(gate, qargs, mapping, connections, compiled_circuit, tot_swaps):
    """
    Gestisce l'applicazione di una singola istruzione sul circuito compilato.
    
    Args:
        gate: Il gate da applicare.
        qargs: Lista dei qubit logici coinvolti.
        mapping: Mappatura dei qubit logici su quelli fisici.
        connections: Lista delle connessioni tra i qubit fisici.
        compiled_circuit: Il circuito compilato su cui aggiungere i gate.
    
    Returns:
        mapping: Mappatura aggiornata dei qubit logici su quelli fisici.
    """
    if len(qargs) == 1:
        # Single-qubit operation
        physical_qubit = mapping[qargs[0]._index]
        # print(f"Applying {gate.name} on physical qubit {physical_qubit}.")
        compiled_circuit.append(gate, [physical_qubit])  # Aggiungi l'operazione al circuito compilato
    
    elif len(qargs) == 2:
        # Two-qubit operation
        logical_q1 = qargs[0]._index
        logical_q2 = qargs[1]._index
        physical_q1 = mapping[logical_q1]
        physical_q2 = mapping[logical_q2]
        
        if isConnected(connections, physical_q1, physical_q2):
            # Se i qubit fisici sono connessi, applica direttamente
            # print(f"Applying {gate.name} on physical qubits {physical_q1} and {physical_q2}.")
            compiled_circuit.append(gate, [physical_q1, physical_q2])  # Aggiungi il gate al circuito compilato
        else:
            # Se i qubit fisici non sono connessi, esegui il routing
            #print(f"Routing required for {gate.name} between logical qubit {logical_q1} (mapped to physical {physical_q1}) and logical qubit {logical_q2} (mapped to physical {physical_q2}).")
            tot_swaps += applyRouting(connections, physical_q1, physical_q2, mapping, compiled_circuit)
            # Aggiorna mapping e ottieni i nuovi qubit fisici          
            routed_q1 = mapping[logical_q1]
            routed_q2 = mapping[logical_q2]
            #print(f"Routed easily following the map chat {routed_q1} {routed_q2}")
            
            # Applica il gate dopo il routing
            #print(f"Applying {gate.name} on updated physical qubits {routed_q1} and {routed_q2}.")
            compiled_circuit.append(gate, [routed_q1, routed_q2])  # Aggiungi il gate al circuito compilato
            
            #print(f"Adesso il mapping è: {mapping}")
    
    else:
        raise ValueError(f"Unsupported gate with {len(qargs)} qubits: {gate.name}")
    
    return mapping, tot_swaps

def bfs(connections, start, end):
    """Trova il percorso minimo tra start e end usando BFS."""
    # La coda per BFS
    queue = deque([[start]])
    visited = set([start])
    
    while queue:
        path = queue.popleft()
        current = path[-1]
        
        if current == end:
            return path  # Restituisce il percorso
       
        for neighbor in connections[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append(new_path)
    
    return []  # Restituisce un percorso vuoto se non c'è connessione

def applyRouting(connections, q1, q2, mapping, compiled_circuit):
    """
    Applica il routing per connettere due qubit fisici non direttamente connessi
    e mantiene il mapping iniziale.
    
    Args:
        connections: Lista delle connessioni tra i qubit fisici.
        q1, q2: Indici dei qubit fisici da connettere.
        mapping: Mappatura dei qubit logici su quelli fisici.
        compiled_circuit: Il circuito compilato su cui aggiungere gli SWAP.
    
    Returns:
        'n_swaps': number of swap operations required
    """
    # Trova il percorso minimo da q1 a q2 
    path = bfs(connections, q1, q2) 
     
    n_swaps = 0
    if not path:
        raise ValueError(f"No path found between qubits {q1} and {q2}.")
    
    # Fai lo swap lungo il percorso
    for i in range(len(path) - 2):   # con -2 (anziché -1) mi sto fermando ad un "vicino del target"
        qubit_j = path[i]
        qubit_k = path[i+1]
        # Aggiungi SWAP tra i qubit lungo il percorso
        # compiled_circuit.swap(path[i], path[i + 1])   # swap ---> native gates
        
        # Implementa lo swap con CX
        compiled_circuit.cx(qubit_j, qubit_k)
        compiled_circuit.cx(qubit_k, qubit_j)
        compiled_circuit.cx(qubit_j, qubit_k)   # swap non è un native gate, ma cx sì
        
        n_swaps += 1
        
        # Aggiorna il mapping (swap indici)
        logic_qbit_j = mapping.index(qubit_j)
        logic_qbit_k = mapping.index(qubit_k)
        mapping[logic_qbit_j], mapping[logic_qbit_k] = mapping[logic_qbit_k], mapping[logic_qbit_j]
        #print(f"Iteration {i}. Update mapping: {mapping}")   

    return n_swaps

def reduceStatevector(statevector, mapping, num_qubits):
    num_states = 2**num_qubits     # Numero totale di stati possibili
    data = [0+0j] * num_states    # Array per memorizzare lo statevector ridotto
    
    # Per ogni stato non nullo dello statevector 
    tot_states = 2**max_qubits
    for index in range(tot_states):
        state = statevector[index]
        if state != 0:      # Considera solo stati non nulli
            # Calcola il nuovo indice ridotto in base al mapping
            new_index = 0
            for logical_qubit in range(num_qubits):
                physical_qubit = mapping[logical_qubit]  # Mappatura da logico a fisico
                if (index >> physical_qubit) & 1:       # Controlla il bit corrispondente
                    new_index += 2**logical_qubit       # Aggiorna l'indice ridotto
            
            # Aggiorna il valore nello statevector ridotto
            data[new_index] += statevector[index]
            #data[new_index] += statevector.data[index] 

    # Restituisce il nuovo statevector con la permutazione applicata
    return Statevector(data, dims=[2] * num_qubits)

# Main
if __name__ == "__main__":
    # Loop sui file QASM
    for name in file_names:
        filename = os.path.join(ABS_PATH, name)
        print(f"Processing file: {filename}")
        
        circuit, translated_circuit, num_qubits = process_qasm_file(filename, native_gates)
        nonCompDepth = circuit.depth()
        nonCompGateCount = sum(circuit.count_ops().values())
        statevector = Statevector(circuit)
        
        print("Original circuit:")
        #print(circuit.draw())
        #plt.show()
        
        #         ###### TRIAL ####
        # # Converti lo Statevector in un array numpy
        # statevector_array = np.array(statevector.data)

        # # Trova gli indici dei valori non nulli
        # nonzero_indices = np.nonzero(statevector_array)[0]

        # # Recupera i valori non nulli
        # nonzero_values = statevector_array[nonzero_indices]

        # # Stampa i risultati
        # print("Indici dei valori non nulli:", nonzero_indices)
        # print("Valori non nulli:", nonzero_values)
        
        # # Converte gli indici in rappresentazioni binarie
        # num_qubitsStatevector = statevector.num_qubits
        # print("Num qubits original stavector: ", num_qubitsStatevector)
        # print("Num qubits original circuit: ", num_qubits)

        # binary_indices = [format(i, f'0{num_qubits}b') for i in nonzero_indices]

        # print("Indici binari:", binary_indices)
                
        # #### and ERROR

        # 3) copy the circuit, add the measurement and simulate with the BasicProvider;
        qc_m = circuit.copy()         # copy the circuit
        qc_m.measure_all()       # add the measurement

        result = backend_basic.run(qc_m).result()
        counts = result.get_counts()
        plot_histogram(counts)
        # Generate all possible binary outcomes for num_qubits
        all_states = [''.join(state) for state in product('01', repeat=num_qubits)]

        # Calculate probabilities, including zero counts
        total_shots = sum(counts.values())
        probabilities = [counts.get(state, 0) / total_shots for state in all_states]
  
        # Salvare la topologia delle connessioni. backend.instructions dizionario 
        # ---> Una lista di liste: in posizione i salvi tutti gli elementi connessi a quel qubit fisico
        connections = allowedConnections(backend)   
      
        # Trivial mapping
        mapping = list(range(max_qubits))  # Logical-to-physical qubit mapping
        
        # Processa il circuito
        compiled_circuit = QuantumCircuit(max_qubits)
        tot_swaps = 0   # Additional swaps required for routing
        print("Processing the circuit...")
        for instr in translated_circuit.data:
            gate = instr.operation  # Access the gate operation
            qargs = instr.qubits  # Access the qubits
            mapping, tot_swaps = apply_gate(gate, qargs, mapping, connections, compiled_circuit, tot_swaps)

        print(f"Final mapping: {mapping}")
        
        # print("Compiled circuit:")
        #print(compiled_circuit.draw())
        # plt.show()
        
        depth = compiled_circuit.depth()
        gateCount = sum(compiled_circuit.count_ops().values())
        
        statevector2 = Statevector(compiled_circuit)
        # print("Stavector: \n", statevector)
        print("Stavector compiled: \n", statevector2)

        #### metodo "probabilities"
        # statevector3 = statevector2.probabilities(mapping[:num_qubits])
        # print("Stavector compiled and reduced (documentation): \n", statevector3)

        statevector3 = reduceStatevector(statevector2, mapping, num_qubits)
        print("Stavector compiled and reduced: \n", statevector3)
        
        #statevector_fidelity = None
        statevector_fidelity = np.abs(statevector.inner(statevector3))
        print(f"Fidelity between non compiled and compiled statevectors: {statevector_fidelity} \n")
        
        # Measure  
        compiled_circuit.measure_all()      # add the measurement

        result2 = backend.run(compiled_circuit).result()
        counts2 = result.get_counts()

        # Calculate probabilities, including zero counts
        num_qubits2 = compiled_circuit.num_qubits
        total_shots2 = sum(counts2.values())
        all_states2 = [''.join(state) for state in product('01', repeat=num_qubits2)]
        probabilities2 = [counts2.get(state, 0) / total_shots2 for state in all_states]
        print("Counts original: ", counts)
        print(f"Counts of the compiled circuit: \n{counts2}")
                
        # To compute fidelity between two statevectors
        prob_fidelity = (np.sum(np.sqrt(probabilities) * np.sqrt(probabilities2)))**2
        print(f"Fidelity between non compiled and compiled probabilities: {prob_fidelity} \n")
       
        # Aggiungi i risultati direttamente al DataFrame
        df.loc[len(df)] = [name, nonCompDepth, nonCompGateCount, depth, gateCount, tot_swaps, prob_fidelity, statevector_fidelity] 

    # Verifica se il file esiste ed è vuoto
    if not os.path.exists(csv_file) or os.path.getsize(csv_file) == 0:
        # Se il file non esiste o è vuoto, scrivi l'header
        df.to_csv(csv_file, mode='w', header=True, index=False)
    else:
        # Se il file esiste ed ha già dei dati, scrivi solo i dati senza l'header
        df.to_csv(csv_file, mode='a', header=False, index=False)

    print("Data successfully saved in the file es02.csv")
        
            
    # Compute total time of execution
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Total time of execution: {total_time:.2f} seconds")
