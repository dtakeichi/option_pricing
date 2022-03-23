import copy
import numpy as np

"""
図3.13
"""
def implicit_FDM_american_put(K, T, S, sig, r, div, N, Nj, dx):
    """
    Parameters
    ----------
    K : float
        行使価格.
    T : float
        満期.
    S : float
        原資産価格.
    sig : float
        ボラティリティ.
    r : float
        リスクフリーレート.
    div : float
        配当率.
    N : float
        タイムステップ数.
    dx : float
        変動幅.

    Returns
    陰的差分法によるアメリカンプット
    使用例 implicit_FDM_american_put(100, 1, 100, 0.2, 0.06, 0.03, 3, 3, 0.2)
    """
    dt = T/N
    nu = r - div - 0.5*sig**2
    edx = exp(dx)
    pu = -0.5*dt*((sig/dx)**2 + nu/dx)  # 式(3.18)
    pd = -0.5*dt*((sig/dx)**2 - nu/dx) 
    pm = 1.0 + dt*(sig/dx)**2 + r*dt
    
    # 満期タイムステップNにおける原資産価格の初期設定
    St = [0 for i in range(2*Nj+1)]
    St[0] = S*exp(Nj*dx)
    for j in range(1,2*Nj+1):
        St[j] = St[j-1]/edx
    St = np.array(St)
    
    # 満期によるオプション価格の初期設定
    C = [[0.0 for i in range(0, 2*Nj+1)] for i in range(0,N+1)]
    C = np.array(C)
    
    # 満期Nのペイオフ 
    C[N] = np.maximum(0, K - St) 
    
    # 境界条件
    lambda_U = 0.0
    lambda_L =  -1 * (St[2*Nj-1] - St[2*Nj])
    
    # pu, pm, pdからなる行列を生成 (3.40)式
    P = [[0 for i in range(0, 2*Nj+1)] for i in range(0,2*Nj+1)]
    P[0][0] = 1
    P[0][1] = -1
    P[2*Nj][2*Nj] = -1
    P[2*Nj][2*Nj-1] = 1
    for i in range(1, 2*Nj):
        P[i][i-1] = pu
        P[i][i] = pm
        P[i][i+1] = pd
    P = np.array(P)
    
    # ラティスのバックステップ
    for i in range(0,N)[::-1]: # i: N-1 -> 0
        # C[i]について連立方程式を解く (3重対角行列方程式)
        # P*C[i] = C[i+1]
        y = copy.copy(C[i+1])
        y[0] = lambda_U
        y[-1] = lambda_L
        C[i] = np.linalg.solve(P,y)
        # 早期行使条件の適用
        C[i] = np.maximum(C[i], K - St)
        # print(C)
    return C[0][Nj]
