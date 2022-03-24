def negative_corr_and_cv_europearn_call(K, T, S, sig, r, div, N, M):
    """
    p123 負の相関法と制御変量法の組み合わせ
    Parameters
    ----------
    K : float  行使価格.
    T : float  満期.
    S : float  原資産価格.
    sig : float ボラティリティ.
    r : float   リスクフリーレート.
    div : float 配当率.
    N : float  タイムステップ数.
    M : float   試行回数.
    
    Returns
    使用例 nc_and_cv_europearn_call(100, 1, 100, 0.2, 0.06, 0.03, 10, 100)
    """
    dt = T/N
    nudt = (r - div - 0.5*sig**2) * dt
    sigsdt = sig*np.sqrt(dt)
    erddt = exp((r - div)*dt)
    beta1 = -1
    
    sum_CT = 0
    sum_CT2 = 0
    
    for j in range(M):
        St1 = S
        St2 = S
        cv1 = 0 #制御変量
        cv2 = 0
        for i in range(N):
            t = i*dt 
            delta1 = Black_Scholes_delta(S=St1, K=K, T=T, sig=sig, r=r, div=div)
            delta2 = Black_Scholes_delta(S=St2, K=K, T=T, sig=sig, r=r, div=div)
            epsilon = np.random.randn()
            Stn1 = St1*exp(nudt + sigsdt*epsilon)
            Stn2 = St2*exp(nudt + sigsdt*-epsilon)
            cv1 = cv1 + delta1*(Stn1 - St1*erddt)
            cv2 = cv2 + delta2*(Stn2 - St2*erddt)
            St1 = Stn1
            St2 = Stn2
        
        CT = 0.5*max(0, St1 - K) + beta1*cv1 + 0.5*max(0, St2 - K) + beta1*cv2 
        sum_CT = sum_CT + CT
        sum_CT2 = sum_CT2 + CT**2
        
    call_value = sum_CT/M * exp(-r*T)
    SD = np.sqrt((sum_CT2 - sum_CT**2/M) *exp(-2*r*T)/(M-1))
    SE = SD/np.sqrt(M)
    print(call_value)
    print(SD)
    print(SE)

    
    
################################
def Black_Scholes_delta(K, T, S, sig, r, div, N=3, Nj=3, dx=0.2):
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

    C1 = C[0][N-1]
    C2 = C[0][N+1]
    S1 = St[N-1]
    S2 = St[N+1]

    return (C1-C2)/(S1-S2)
