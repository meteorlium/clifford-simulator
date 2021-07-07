from collections import Counter
from tableau import Tableau
import numpy as np


def measure_tableau(t: Tableau, measure_times=1000):

    n_bit = t.n
    t_copy = Tableau(n_bit)

    results = []
    for _ in range(measure_times):
        t_copy.x = t.x.copy()
        t_copy.z = t.z.copy()
        t_copy.r = t.r.copy()

        measure = [t_copy.measure(i_bit) for i_bit in range(n_bit)]
        result = '|'
        for i_bit in range(n_bit-1, 0-1, -1):
            result += str(measure[i_bit])
        result += '>'
        results.append(result)
    return results


def print_measure_results(results, measure_times=1000):
    print('measure results:')
    # print('format: |q[n] q[n-1] ... q[0]>')
    counter = Counter(results)
    for i in counter:
        print(f'{i}: {counter[i] / measure_times}')


def print_state_string(states):
    print(f'{len(states)} nonzero basis states:')
    for state in states:
        print(state)


def get_state_vector(states):
    def get_single_state(state):
        phase_string = state[:2]
        if phase_string == ' +':
            phase = 1.
        elif phase_string == '+i':
            phase = 1j
        elif phase_string == ' -':
            phase = -1.
        elif phase_string == '-i':
            phase = -1j
        vector_string = state[3:-1]
        vector = np.array([1.])
        for qubit in vector_string:
            if qubit == '1':
                vector = np.kron(np.array([0., 1.]), vector)
            else:
                vector = np.kron(np.array([1., 0.]), vector)
        vector = vector * phase
        return vector

    n_state = len(states)
    n_qubit = len(states[0]) - 4  # e.g. n_qubit of '+i|00>' is 2
    vector = np.zeros(2 ** n_qubit)
    for state in states:
        vector = vector + get_single_state(state)
    vector /= np.sqrt(n_state)
    return vector.reshape([-1, 1])
