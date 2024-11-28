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
from collections import deque

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
max_qubits = backend.configuration().num_qubits    # Ottieni il numero di qubit disponibili dal backend


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
def process_qasm_file(qasm_file, native_gates):   # il numero massimo di qubit è fissato dal backend (16 per fakeguadalupe)
    """Leggi un file QASM, verifica i gate e traduci quelli non nativi."""
    # Carica il circuito dal file QASM
    circuit = QuantumCircuit.from_qasm_file(qasm_file)
    num_qubits = circuit.num_qubits      # è già un parametro di dominio pubblico, che posso chiamare anche fuori dalla funzione?
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
    
    return translated_circuit, num_qubits

# 
def allowedConnections(backend):
    """Crea una lista di liste per rappresentare le connessioni tra i qubit fisici."""
    coupling_map = backend.configuration().coupling_map  # Ottieni la mappa delle connessioni
    num_qubits = backend.configuration().num_qubits      # Numero di qubit nel backend

    # Inizializza la lista delle connessioni con set per rimuovere i duplicati automaticamente
    connections = [set() for _ in range(num_qubits)]

    for edge in coupling_map:
        q1, q2 = edge
        connections[q1].add(q2)
        connections[q2].add(q1)  # Connessione bidirezionale

    # Converti i set in liste per la compatibilità con il resto del codice
    return [list(conn) for conn in connections]

def isConnected(connections, q1, q2):
    """Controlla se due qubit fisici sono direttamente connessi."""
    return q2 in connections[q1]    # 1 se sono connessi

def apply_gate(gate, qargs, mapping, connections, compiled_circuit):
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
        physical_qubit = mapping[qargs[0].index]
        print(f"Applying {gate.name} on physical qubit {physical_qubit}.")
        compiled_circuit.append(gate, [physical_qubit])  # Aggiungi l'operazione al circuito compilato
    
    elif len(qargs) == 2:
        # Two-qubit operation
        physical_q1 = mapping[qargs[0].index]
        physical_q2 = mapping[qargs[1].index]
        
        if isConnected(connections, physical_q1, physical_q2):
            # Se i qubit fisici sono connessi, applica direttamente
            print(f"Applying {gate.name} on physical qubits {physical_q1} and {physical_q2}.")
            compiled_circuit.append(gate, [physical_q1, physical_q2])  # Aggiungi il gate al circuito compilato
        else:
            # Se i qubit fisici non sono connessi, esegui il routing
            print(f"Routing required for {gate.name} between {physical_q1} and {physical_q2}.")
            routing_result = applyRouting(connections, physical_q1, physical_q2, mapping, compiled_circuit)
            # Aggiorna mapping e ottieni i nuovi qubit fisici
            mapping = routing_result['updated_mapping']    # l'idea è di lasciare il trivial mapping, quindi non varia
            routed_q1 = routing_result['q1']
            routed_q2 = routing_result['q2']
            
            # Applica il gate dopo il routing
            print(f"Applying {gate.name} on updated physical qubits {routed_q1} and {routed_q2}.")
            compiled_circuit.append(gate, [routed_q1, routed_q2])  # Aggiungi il gate al circuito compilato
    
    else:
        raise ValueError(f"Unsupported gate with {len(qargs)} qubits: {gate.name}")
    
    return mapping

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
        dict: Contiene:
            - 'q1', 'q2': I nuovi indici fisici dei qubit dopo il routing.
            - 'updated_mapping': La mappatura aggiornata dei qubit logici.
    """
    # Trova il percorso minimo da q1 a q2
    path = bfs(connections, q1, q2)
    
    if not path:
        raise ValueError(f"No path found between qubits {q1} and {q2}.")
    
    # Fai lo swap lungo il percorso
    for i in range(len(path) - 1):
        # Aggiungi SWAP tra i qubit lungo il percorso
        compiled_circuit.swap(path[i], path[i + 1])
        
        # Aggiorna il mapping
        for logical_qubit, physical_qubit in enumerate(mapping):
            if physical_qubit == path[i]:
                mapping[logical_qubit] = path[i + 1]
            elif physical_qubit == path[i + 1]:
                mapping[logical_qubit] = path[i]
    
    # Ora q1 e q2 sono connessi, quindi restituiamo la mappatura aggiornata
    return {'q1': q1, 'q2': q2, 'updated_mapping': mapping}    

# Memorizza i risultati
def save_results(file_name, original_circuit, compiled_circuit):
    # Calcolare la profondità e il numero di gate nel circuito originale
    original_depth = original_circuit.depth()
    original_gate_count = len(original_circuit.data)
    
    # Calcolare la profondità e il numero di gate nel circuito compilato
    compiled_depth = compiled_circuit.depth()
    compiled_gate_count = len(compiled_circuit.data)
    
    # Calcolare la fedeltà tra il circuito originale e compilato
    original_state = Statevector.from_instruction(original_circuit)
    compiled_state = Statevector.from_instruction(compiled_circuit)
    statevector_fidelity_value = state_fidelity(original_state, compiled_state)
    
    # Includere la fedeltà delle probabilità
    # Se vuoi includere anche la fedeltà basata sulle probabilità, puoi simulare i circuiti e calcolare la fidelità sulle probabilità
    original_probabilities = original_state.probabilities_dict()
    compiled_probabilities = compiled_state.probabilities_dict()
    
    probability_fidelity_value = sum(min(original_probabilities.get(key, 0), compiled_probabilities.get(key, 0)) 
                                     for key in set(original_probabilities) | set(compiled_probabilities))
    
    # Aggiungere i risultati alla lista
    results.append({
        "Name of the circuit": file_name,
        "Original Depth": original_depth,
        "Original Gate Count": original_gate_count,
        "Compiled Depth": compiled_depth,
        "Compiled Gate Count": compiled_gate_count,
        "Statevector Fidelity": statevector_fidelity_value,
        "Probability Fidelity": probability_fidelity_value
    })



# Main
# Loop sui file QASM
for name in file_names:
    filename = os.path.join(ABS_PATH, name)
    print(f"Processing file: {filename}")
    
    translated_circuit, num_qubits = process_qasm_file(filename, native_gates)
    print(translated_circuit)
    
    # Salvare la topologia delle connessioni. backend.instructions dizionario 
    # ---> Una lista di liste: in posizione i salvi tutti gli elementi connessi a quel qubit fisico
    connections = allowedConnections(backend)     # il grafo è uguale per ogni singola operazione? probabilmente i single qubit erano possibili in ogni posizione. A noi serve capire
    
    # Trivial mapping (e lo manteniamo fino alla fine)
    mapping = list(range(num_qubits))  # Logical-to-physical qubit mapping
    
    # Processa il circuito
    compiled_circuit = QuantumCircuit(max_qubits)
    
    print("Processing the circuit...")
    for instr in translated_circuit.data:
        gate = instr[0]
        qargs = instr[1]
        mapping = apply_gate(gate, qargs, mapping, connections, compiled_circuit)

    print("Compiled circuit:")
    print(compiled_circuit)

    # Calcola e salva i risultati
    save_results(name, translated_circuit, compiled_circuit)

# Salvare i risultati in un file CSV
df = pd.DataFrame(results)
df.to_csv(csv_file, index=False)
print(f"Results saved to {csv_file}")
        # se è un single-qubit operation, applicala
    # Prima di applicare l'istruzione "compilata", 
    # isConnected()  che ti dice se è un'operazione lecita (sono connessi i due qubit di una cx): restituisce 1 (True)
     #applica istruzione sui qubit fisici 
        # cx(0,1)   ---> cx(map[0],map[1])  perché ora il qubit logico 0 è mappato nel qubit fisico map[i]

    # Altrimenti routing: SWAP(finché non avvicini i due qubit) e scambiare l'array
    # Una volta applicato il gate, riapplica al contrario gli swap in modo da ripristinare il "trivial mapping"
    
    
    
    # Fidelity check: abbiamo fatto un buon lavoro?
    # Validate it using the provided qasm file, using the fidelity compared to the not compiled circuit results.
    
    # Basta misurare i primi (o ultimi?) num_qubits  [trivial mapping]
    # Fidelity check
        
# Compute total time of execution
end_time = time.time()
total_time = end_time - start_time
print(f"Total time of execution: {total_time:.2f} seconds")


    


"""
    # Memorizzare i risultati
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
       



    
