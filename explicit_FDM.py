import numpy as np
from math import exp
from math import sqrt

def explicit_FDM_european_call(K, T, S, sig, r, div, N, Nj, dx):
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
    陽的差分法によるヨーロピアンコール
    使用例 explicit_FDM_european_call(100, 1, 100, 0.2, 0.06, 0.03, 3, 3, 0.2)
    """
    dt = T/N
    nu = r - div - 0.5*sig**2
    edx = exp(dx)
    pu = 0.5*dt*((sig/dx)**2 + nu/dx)  # 式(3.18)
    pd = 0.5*dt*((sig/dx)**2 - nu/dx) 
    pm = 1.0 - dt*(sig/dx)**2 - r*dt
    
    # 満期タイムステップNにおける原資産価格の初期設定
    St = [0 for i in range(2*Nj+1)]
    St[0] = S*exp(-Nj*dx)
    for j in range(1,2*Nj+1):
        St[j] = St[j-1]*edx
    
    # 満期によるオプション価格の初期設定
    C = [[0 for i in range(0, 2*Nj+1)] for i in range(0,N+1)]
    # 満期Nのペイオフ
    for j in range(0,2*Nj+1): 
        C[N][j] = max(0, St[j] - K)
    
    # ラティスのバックステップ
    for i in range(0,N)[::-1]: # i: N-1 -> 0
        for j in range(1,2*Nj):
            C[i][j] = pu*C[i+1][j+1] + pm*C[i+1][j] + pd*C[i+1][j-1]   
        # 境界条件
        C[i][0] = C[i][1]
        C[i][2*Nj] = C[i][2*Nj-1] + (St[2*Nj] - St[2*Nj-1])

    return C[0][N]
    
 

