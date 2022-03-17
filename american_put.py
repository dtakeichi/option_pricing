from math import exp

def american_put(K, T, S, r, N, u, d):
    """
    Parameters
    ----------
    K : float
        行使価格.
    T : float
        満期.
    S : float
        原資産価格.
    r : TYPE
        リスクフリーレート.
    N : float
        タイムステップ数.
    u : float
        上昇率.
    d : float
        下降率.

    Returns
    アメリカンプット
    """
    dt = T/N
    p = (exp(r*dt)-d) / (u-d) #リスク中立の上昇確率
    disc = exp(-r*dt) # 割引率
    
    # 満期タイムステップNにおける原資産価格の初期設定
    St = [0 for i in range(N+1)]
    St[0] = S*(d**N)
    
    for j in range(1,N+1):
        St[j] = St[j-1]*u/d
        
    # 満期によるオプション価格の初期設定
    C = [0 for i in range(0,N+1)]
    
    for j in range(0,N+1):
        C[j] = max(0, K - St[j]) #プット
        
    # 早期行使条件を適用しながらツリーのステップバック
    for i in range(0,N)[::-1]: # i: N-1 -> 0
        for j in range(0,i+1): # j: 0 -> i
            C[j] = disc*(p*C[j+1] + (1-p)*C[j])
            St[j] = St[j]/d
            C[j] = max(C[j], K - St[j])

    return C[0]
