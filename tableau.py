import numpy as np
import random


class Tableau(object):
    """simulation of Clifford circuits

    Reference:
        "Improved Simulation of Stabilizer Circuits"
        Scott Aaronson and Daniel Gottesman
        https://arxiv.org/abs/quant-ph/0406196
    """

    def __init__(self, num_qubits):
        self.n = num_qubits
        tableau = np.eye(2 * num_qubits + 1, dtype=np.int8)
        # row 2n+1 is scratch space, used for support computation
        # x_ij, z_ij in {0, 1}
        self.x = tableau[:, :self.n]
        self.z = tableau[:, self.n: self.n*2]
        # r_i in {0, 1, 2, 3}. means phase = 1j^r.
        self.r = np.zeros((2 * num_qubits + 1, 1), dtype=np.int8)

    # clifford gates

    def cx(self, a, b):
        """CNOT from control a to target b
        """
        r_temp = self.x[:, a] & self.z[:, b] & (
            self.x[:, b] ^ self.z[:, a] ^ True)
        self.r[:, 0] = (self.r[:, 0] + r_temp * 2) % 4
        self.x[:, b] ^= self.x[:, a]
        self.z[:, a] ^= self.z[:, b]

    def h(self, a):
        """Hadamard on qubit a
        """
        r_temp = self.x[:, a] & self.z[:, a]
        self.r[:, 0] = (self.r[:, 0] + r_temp * 2) % 4
        temp = self.x[:, a].copy()
        self.x[:, a] = self.z[:, a].copy()
        self.z[:, a] = temp

    def s(self, a):
        """Phase on qubit a
        """
        r_temp = self.x[:, a] & self.z[:, a]
        self.r[:, 0] = (self.r[:, 0] + r_temp * 2) % 4
        self.z[:, a] ^= self.x[:, a]

    # row operation

    def _row_mul(self, a, b):
        """
        row_a表示矩阵A
        row_b表示矩阵B
        计算B*A，保存在row_a
        """
        def row_mul_r(x1, z1, x2, z2):
            """
            论文中的g函数
            计算矩阵乘积后是否有符号变化，即i或-i
            """
            if not x1 and not z1:  # I
                return 0
            if x1 and z1:  # Y
                return z2 - x2
            if x1 and not z1:  # X
                return z2 * (2 * x2 - 1)
            if not x1 and z1:  # Z
                return x2 * (1 - 2 * z2)
        r_change = 0
        for jj in range(self.n):
            r_change += row_mul_r(self.x[b, jj], self.z[b, jj],
                                  self.x[a, jj], self.z[a, jj])

        self.r[a] = (self.r[a] + self.r[b] + r_change) % 4
        self.x[a, :self.n] ^= self.x[b, :self.n]
        self.z[a, :self.n] ^= self.z[b, :self.n]

    def _row_swap(self, a, b):
        """swap row a and b
        """
        temp = self.x[a].copy()
        self.x[a] = self.x[b].copy()
        self.x[b] = temp

        temp = self.z[a].copy()
        self.z[a] = self.z[b].copy()
        self.z[b] = temp

        temp = self.r[a].copy()
        self.r[a] = self.r[b].copy()
        self.r[b] = temp

    def _clear_scratch_space(self):
        """clear the scratch space, row 2n"""
        self.x[2 * self.n] = 0
        self.z[2 * self.n] = 0
        self.r[2 * self.n] = 0

    # measure

    def measure(self, a):
        """measurement of qubit a in standard basis
        """
        for p in range(self.n, 2 * self.n):
            if self.x[p, a] == 1:
                # case 1: 0或1都有可能，随机决定并改变tableau
                return self._measure_random(a, p)
        # case 2: 有确定结果
        return self._measure_determinate(a)

    def _measure_random(self, a, p):
        n = self.n
        for ii in range(2 * n):
            if (ii != p) and self.x[ii, a]:
                self._row_mul(ii, p)
        self.x[p - n] = self.x[p]
        self.z[p - n] = self.z[p]
        self.r[p - n] = self.r[p]
        self.x[p] = 0
        self.z[p] = 0
        result = random.randint(0, 1)
        self.r[p] = result * 2
        self.z[p, a] = 1
        return result

    def _measure_determinate(self, a):
        n = self.n
        self._clear_scratch_space()
        for ii in range(n):
            if self.x[ii, a]:
                self._row_mul(2 * n, ii + n)
        return self.r[2 * n, 0] // 2

    # compute qubit state

    def gaussian(self):
        """chp: gaussian

        将tableau通过行变换进行格式化，类似矩阵的高斯消去。

        目标：
        将x和z的下半部分，即stab set的部分，变成如下格式：
        分成上下两部分，上半部分只有XY和I。下半部分只有Z和I。两部分都变成上三角形。
        消去的方法与高斯消去法类似，但不适用加减而是用乘法。
        原理是，若有矩阵U和V以及向量psi，满足U*psi=psi，V*psi=psi，则有U*V*psi=psi。
        对x和z的下半部分进行格式化的同时对相应的上半部分变化。

        用处：
        是计算量子态(compute_ket函数)的前置步骤。
        """

        n = self.n
        i = n

        for j in range(n):
            for k in range(i, 2 * n):
                if self.x[k, j] == 1:  # Find a generator containing X in jth column
                    self._row_swap(i, k)
                    self._row_swap(i - n, k - n)
                    for k2 in range(i + 1, 2 * n):
                        if self.x[k2, j] == 1:
                            self._row_mul(k2, i)  # Gaussian elimination step
                            self._row_mul(i - n, k2 - n)
                    i += 1
                    break
        g = i - n

        for j in range(n):
            for k in range(i, 2 * n):  # Find a generator containing Z in jth column
                if self.z[k, j] == 1:
                    self._row_swap(i, k)
                    self._row_swap(i - n, k - n)
                    for k2 in range(i + 1, 2 * n):
                        if self.z[k2, j] == 1:
                            self._row_mul(k2, i)  # Gaussian elimination step
                            self._row_mul(i - n, k2 - n)
                    i += 1
                    break
        return g

    def compute_ket(self):
        """chp: print ket
        """
        g = self.gaussian()

        """
        chp.seed
        Finds a Pauli operator P such that the basis state P | 0...0 > occurs with nonzero amplitude in the tableau, and writes P to the scratch space.
        For this to work, Gaussian elimination must already have been performed on the tableau.
        找到矩阵P，使得P*|0...0>=psi为量子态中的一个。将P保存至第2n行。

        return states: list of string.
        e.g. [' +|00>', ' +|11>']
        """
        self._clear_scratch_space()  # wipe the scratch space clean
        n = self.n
        for i in range(2 * n - 1, n + g - 1, -1):
            # gaussian以后，n+g-1到2n-1行只有Z和I，没有XY
            f = self.r[i]
            for j in range(n - 1, 0 - 1, -1):
                if self.z[i, j] == 1:
                    min_z = j
                    if self.x[2 * n, j] == 1:
                        f = (f + 2) % 4
            if f == 2:
                # make the seed consistent with the ith equation
                self.x[2 * n, min_z] = not self.x[2 * n, min_z]

        states = []
        states.append(self._get_basis_state())
        for t in range(2 ** g - 1):
            t2 = t ^ (t + 1)
            for i in range(g):
                if t2 & (2 ** i):
                    self._row_mul(2 * n, n + i)
            states.append(self._get_basis_state())
        return states

    def _get_basis_state(self):
        """chp: print basis state
        Prints the result of applying the Pauli operator in the "scratch space" of q to |0...0>
        打印compute-ket函数中保存在第2n行的Pauli矩阵决定的向量。

        return state: str, e.g. +i|001>
        """
        state = ''
        n = self.n
        i_exp = self.r[2 * n]
        for j in range(n):
            if self.x[2 * n, j] == 1 and self.z[2 * n, j] == 1:
                i_exp = i_exp + 1
        i_exp = i_exp % 4
        if i_exp == 0:
            state += ' +'
        elif i_exp == 1:
            state += '+i'
        elif i_exp == 2:
            state += ' -'
        elif i_exp == 3:
            state += '-i'
        state += '|'

        for j in range(n-1, 0-1, -1):  # 倒序
            if self.x[2 * n, j] == 1:
                state += '1'
            else:
                state += '0'
        state += '>'
        return state

    # show

    def show_tableau(self):
        print('tableau:')
        for ii in range(self.n):
            for jj in range(self.n):
                print(self.x[ii, jj], end=' ')
            print('|', end=' ')
            for jj in range(self.n):
                print(self.z[ii, jj], end=' ')
            print('|', end=' ')
            print(self.r[ii, 0])
        print('-' * self.n * 2, end='')
        print('|-', end='')
        print('-' * self.n * 2, end='')
        print('|-', end='')
        print('-')
        for ii in range(self.n, self.n * 2):
            for jj in range(self.n):
                print(self.x[ii, jj], end=' ')
            print('|', end=' ')
            for jj in range(self.n):
                print(self.z[ii, jj], end=' ')
            print('|', end=' ')
            print(self.r[ii, 0])

    def show_stab(self):
        """show stab(psi)
        """
        print('stab set:')
        for ii in range(self.n, self.n * 2):
            if self.r[ii, 0] == 0:
                print(' +', end='')
            elif self.r[ii, 0] == 1:
                print('+i', end='')
            elif self.r[ii, 0] == 2:
                print(' -', end='')
            elif self.r[ii, 0] == 3:
                print('-i', end='')

            for jj in range(self.n):
                x = self.x[ii, jj]
                z = self.z[ii, jj]

                """encode Pauli matrix using x and z
                00: I
                01: X
                11: Y
                10: Z
                """

                if x == 0 and z == 0:
                    print('I', end='')
                elif x == 0 and z == 1:
                    print('Z', end='')
                elif x == 1 and z == 1:
                    print('Y', end='')
                elif x == 1 and z == 0:
                    print('X', end='')
            print()
