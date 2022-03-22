import numpy as np
from math import exp
from math import sqrt

def trinomial_tree_european_call(K, T, S, sig, r, div, N, dx):
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
    3項ツリーによるヨーロピアンコール
    """
    dt = T/N
    nu = r - div - 0.5*sig**2
    edx = exp(dx)
    pu = 0.5*((sig**2*dt + nu**2*dt**2) / dx**2 + nu*dt/dx) #上昇確率
    pd = 0.5*((sig**2*dt + nu**2*dt**2) / dx**2 - nu*dt/dx) #下降確率
    pm = 1.0 - (sig**2*dt + nu**2*dt**2) / dx**2 #変化しない確率

    disc = exp(-r*dt) # 割引率
    
    # 満期タイムステップNにおける原資産価格の初期設定
    St = [0 for i in range(2*N+1)]
    St[0] = S*exp(-N*dx)
    for j in range(1,2*N+1):
        St[j] = St[j-1]*edx
    
    # 満期によるオプション価格の初期設定
    C = [[0 for i in range(0, 2*N+1)] for i in range(0,N+1)]

    for j in range(0,2*N+1):
        C[N][j] = max(0, St[j] - K)
    
    # 時間のバックステップ
    for i in range(0,N)[::-1]: # i: N-1 -> 0
        # N=3 のとき, N=1なら3, N=2なら5なので N=iなら2i+1個ある
        for j in range(0,2*i+1): # j: 0 -> i
            C[i][j] = disc*(pu*C[i+1][j+2] + pm*C[i+1][j+1] + pd*C[i+1][j])

    return C[0][0]
