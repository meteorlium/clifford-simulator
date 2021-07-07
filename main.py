import util
import numpy as np
from tableau import Tableau


def circuit():
    """
    测试线路
    https://algassert.com/quirk#circuit={%22cols%22:[[%22H%22,%22H%22,1,1,%22H%22],[%22%E2%80%A2%22,1,%22X%22],[%22%E2%80%A2%22,1,1,%22X%22],[%22%E2%80%A2%22,1,1,1,1,%22X%22],[%22Z^%C2%BD%22,%22Z^%C2%BD%22,%22Z^%C2%BD%22,1,1,%22Z^%C2%BD%22],[1,1,1,%22X%22,1,%22%E2%80%A2%22]]}
    """
    n_bit = 6
    t = Tableau(n_bit)
    t.h(0)
    t.h(1)
    t.h(4)
    t.cx(0, 2)
    t.cx(0, 3)
    t.cx(0, 5)
    t.s(0)
    t.s(1)
    t.s(2)
    t.s(5)
    t.cx(5, 3)
    return t

# def circuit():
#     """测试线路：ghz without measurement
#     https://algassert.com/quirk#circuit={%22cols%22:[[%22H%22,%22H%22],[%22%E2%80%A2%22,1,%22X%22],[1,%22%E2%80%A2%22,%22X%22],[%22Z^%C2%BD%22,%22Z^%C2%BD%22,%22Z^%C2%BD%22],[%22H%22,%22H%22]]}
#     """
#     n_bit = 3
#     t = Tableau(n_bit)
#     t.h(0)
#     t.h(1)
#     t.cx(0, 2)
#     t.cx(1, 2)
#     t.s(0)
#     t.s(1)
#     t.s(2)
#     t.h(0)
#     t.h(1)
#     return t

# def circuit():
#     n_bit = 2
#     t = Tableau(n_bit)
#     t.h(0)
#     t.cx(0, 1)
#     t.s(1)
#     return t


if __name__ == '__main__':

    # 构造线路t

    t = circuit()

    # 输出tableau表格和稳定算子

    t.gaussian()
    print()
    t.show_tableau()
    print()
    t.show_stab()

    # 计算量子态的ket形式

    states = t.compute_ket()

    print()
    util.print_state_string(states)

    print()
    vector = util.get_state_vector(states)
    print('state vector:')
    print(vector)

    density_matrix = np.matmul(vector, vector.T)
    print()
    print('density matrix:')
    print(density_matrix)

    # 测量
    count_time = True
    if count_time:
        import time
        starttime = time.time()
    measure_times = 10000
    measure_results = util.measure_tableau(t, measure_times)
    print()
    util.print_measure_results(measure_results, measure_times)
    if count_time:
        endtime = time.time()
        print(f'measure time: {endtime - starttime} seconds')
