import networkx as nx
import heapq


def calculate_makespan(tasks, dependencies, num_computers, configure_node_wait_queue_time=60/3600):
    """
    Calculate the makespan of a DAG given task times and dependencies.

    Args:
        tasks (dict): Task times (e.g., {'A': 3, 'B': 2}).
        dependencies (list of tuples): Edges representing dependencies (e.g., [('A', 'B'), ('B', 'C')]).
        num_computers (int): Number of available computers.
        configure_node_wait_queue_time (float): How much time does it take to have the node ready for the next
        job; by default 60 seconds. It must be converted to hours

    Returns:
        int: Total completion time (makespan).
        dict: Task schedule with start and end times.
    """
    # Step 1: Build the DAG
    dag = nx.DiGraph()
    dag.add_nodes_from(tasks.keys())
    dag.add_edges_from(dependencies)

    # Step 2: Topological Sort
    try:
        top_order = list(nx.topological_sort(dag))
    except nx.NetworkXUnfeasible:
        raise ValueError("The graph contains a cycle!")

    # Step 3: Schedule Tasks on Computers
    available_computers = [(0, i) for i in range(num_computers)]  # (next_available_time, computer_id)
    heapq.heapify(available_computers)

    task_schedule = {}
    end_times = {node: 0 for node in dag.nodes}  # End times for all tasks
    current_time = 0

    for node in top_order:
        # Determine earliest start time based on dependencies
        predecessors = list(dag.predecessors(node))
        earliest_start = max([end_times[p] for p in predecessors], default=0)

        # Wait for the earliest available computer
        next_available_time, computer_id = heapq.heappop(available_computers)
        start_time = max(earliest_start, next_available_time)
        # It is also added configure_node_wait_queue_time to simulate
        # the latency of the node when releasing and catching the new job
        end_time = start_time + tasks[node] + configure_node_wait_queue_time

        # Schedule the task
        task_schedule[node] = {
            'computer': computer_id,
            'start': start_time,
            'end': end_time
        }

        # Update end time and computer availability
        end_times[node] = end_time
        heapq.heappush(available_computers, (end_time, computer_id))
        current_time = max(current_time, end_time)

    return current_time, task_schedule


def build_dag_from_fep_simulation(num_ligands, performance_ligand, performance_complex,
                                  equi_sim_time_ligand, fep_sim_time_ligand,
                                  equi_sim_time_complex, fep_sim_time_complex):
    """
    Build a DAG for multiple ligands, each with its own simulation tasks.

    Args:
        num_ligands (int): Number of ligands.
        performance_ligand (float): Performance (ns/day) for ligand tasks.
        performance_complex (float): Performance (ns/day) for complex tasks.
        equi_sim_time_ligand (list of float): Equilibration times (ns) for ligand.
        fep_sim_time_ligand (list of of list of float): FEP times (ns) for ligand.
        equi_sim_time_complex (list of float): Equilibration times (ns) for complex.
        fep_sim_time_complex (list of list of float): FEP times (ns) for complex.

    Returns:
        dict: Tasks with their durations (in hours).
        list: Dependencies between tasks.
    """
    tasks = {}
    dependencies = []
    # Convert performance from ns/day to ns/hour
    performance_complex = performance_complex / 24
    performance_ligand = performance_ligand / 24

    for ligand_idx in range(1, num_ligands + 1):
        ligand_prefix = f"lig{ligand_idx}"

        # Process ligand equilibration
        for i, sim_time in enumerate(equi_sim_time_ligand, start=1):
            task_name = f"{ligand_prefix}_equi_lig{i}"
            task_duration = sim_time / performance_ligand  # Convert ns to days
            tasks[task_name] = task_duration
            if i > 1:
                dependencies.append((f"{ligand_prefix}_equi_lig{i-1}", task_name))

        # Process ligand FEP
        for i, sim_times in enumerate(fep_sim_time_ligand, start=1):
            for j, sim_time in enumerate(sim_times, start=1):
                task_name = f"{ligand_prefix}_fep_lig{i}_step{j}"
                task_duration = sim_time / performance_ligand  # Convert ns to days
                tasks[task_name] = task_duration
                if j == 1:
                    dependencies.append((f"{ligand_prefix}_equi_lig{len(equi_sim_time_ligand)}", task_name))
                else:
                    dependencies.append((f"{ligand_prefix}_fep_lig{i}_step{j-1}", task_name))

        # Process complex equilibration
        for i, sim_time in enumerate(equi_sim_time_complex, start=1):
            task_name = f"{ligand_prefix}_equi_comp{i}"
            task_duration = sim_time / performance_complex  # Convert ns to hours
            tasks[task_name] = task_duration
            if i > 1:
                dependencies.append((f"{ligand_prefix}_equi_comp{i-1}", task_name))

        # Process complex FEP
        for i, sim_time in enumerate(fep_sim_time_complex, start=1):
            for j, sim_time in enumerate(sim_times, start=1):
                task_name = f"{ligand_prefix}_fep_comp{i}_step{j}"

                task_duration = sim_time / performance_complex  # Convert ns to hours
                tasks[task_name] = task_duration
                if j == 1:
                    dependencies.append((f"{ligand_prefix}_equi_comp{len(equi_sim_time_complex)}", task_name))
                else:
                    dependencies.append((f"{ligand_prefix}_fep_comp{i}_step{j-1}", task_name))

    return tasks, dependencies


def build_dag_from_mmgbsa_simulation(num_ligands, performance_complex,
                                     equi_sim_time_complex, mmgbsa_sim_time_complex):
    """
    Build a DAG for multiple ligands, each with its own simulation tasks.

    Args:
        num_ligands (int): Number of ligands.
        performance_complex (float): Performance (ns/day) for complex tasks.
        fep_sim_time_ligand (list of float): FEP times (ns) for ligand.
        equi_sim_time_complex (list of float): Equilibration times (ns) for complex.
        mmgbsa_sim_time_complex (list of float): FEP times (ns) for complex.

    Returns:
        dict: Tasks with their durations (in hours).
        list: Dependencies between tasks.
    """
    tasks = {}
    dependencies = []
    # Convert performance from ns/day to ns/hour
    performance_complex = performance_complex / 24

    for ligand_idx in range(1, num_ligands + 1):
        ligand_prefix = f"lig{ligand_idx}"

        # Process complex equilibration
        for i, sim_time in enumerate(equi_sim_time_complex, start=1):
            task_name = f"{ligand_prefix}_equi_comp{i}"
            task_duration = sim_time / performance_complex  # Convert ns to hours
            tasks[task_name] = task_duration
            if i > 1:
                dependencies.append((f"{ligand_prefix}_equi_comp{i-1}", task_name))

        # Process complex FEP
        for i, sim_time in enumerate(mmgbsa_sim_time_complex, start=1):
            task_name = f"{ligand_prefix}_sample_comp{i}"
            task_duration = sim_time / performance_complex  # Convert ns to hours
            tasks[task_name] = task_duration
            dependencies.append((f"{ligand_prefix}_equi_comp{len(equi_sim_time_complex)}", task_name))

    return tasks, dependencies


# Function to calculate makespan for a given number of ligands and computers for fep
def compute_makespan_fep(ligs, comps,
                         performance_ligand, performance_complex,
                         equi_sim_time_ligand,
                         fep_sim_time_ligand,
                         equi_sim_time_complex,
                         fep_sim_time_complex,
                         configure_node_wait_queue_time):
    tasks, dependencies = build_dag_from_fep_simulation(
        ligs,
        performance_ligand, performance_complex,
        equi_sim_time_ligand, fep_sim_time_ligand,
        equi_sim_time_complex, fep_sim_time_complex)
    makespan, _ = calculate_makespan(tasks, dependencies, comps, configure_node_wait_queue_time=configure_node_wait_queue_time)
    return ligs, comps, makespan


# Function to calculate makespan for a given number of ligands and computers for mmgbsa
def compute_makespan_mmgbsa(ligs, comps,
                            performance_complex,
                            equi_sim_time_complex,
                            mmgbsa_sim_time_complex,
                            configure_node_wait_queue_time):
    tasks, dependencies = build_dag_from_mmgbsa_simulation(
        ligs,
        performance_complex,
        equi_sim_time_complex,
        mmgbsa_sim_time_complex)
    makespan, _ = calculate_makespan(tasks, dependencies, comps, configure_node_wait_queue_time=configure_node_wait_queue_time)
    return ligs, comps, makespan