import numpy as np
from math import exp
from math import sqrt
# p56, p58

def american_spread_call(K, T, S1, S2, sig1, sig2, div1, div2, rho, r, N):
    """
    Parameters
    ----------
    K : float
        行使価格.
    T : float
        満期.
    S1 : float
        原資産2価格.
    S2 : float
        原資産1価格.
    sig1 : float
        原資産1標準偏差.
    sig2 : float
        原資産2標準偏差.
    div1 : float
        原資産1の連続的配当
    div2 : float
        原資産2の連続的配当
    rho : float
        2つのブラウン運動の相関係数
    r : float
        リスクフリーレート.

    Returns
    アメリカン・スプレッド・コールオプション

    使用例
    american_spread_call(1, 1, 100, 100, 0.2, 0.3, 0.03, 0.04, 0.5, 0.06, 3) 
    """
    dt = T/N
    nu1 = r - div1 - 0.5*sig1**2 # リスク中立確率のもとでのlogS1のドリフト(配当分低下)
    nu2 = r - div2 - 0.5*sig2**2 # リスク中立確率のもとでのlogS1のドリフト(配当分低下)
    dx1 = sig1*sqrt(dt) # 増分
    dx2 = sig2*sqrt(dt) 
    disc = exp(-r*dt) # 割引率

    #リスク中立の上昇確率の計算 (1,2)でuu, ud, du, ddの4通り
    puu = (dx1*dx2 + (dx2*nu1 + dx1*nu2 + rho*sig1*sig2)*dt) / (4*dx1*dx2) * disc
    pud = (dx1*dx2 + (dx2*nu1 - dx1*nu2 - rho*sig1*sig2)*dt) / (4*dx1*dx2) * disc
    pdu = (dx1*dx2 + (-dx2*nu1 + dx1*nu2 - rho*sig1*sig2)*dt) / (4*dx1*dx2) * disc
    pdd = (dx1*dx2 + (-dx2*nu1 - dx1*nu2 + rho*sig1*sig2)*dt) / (4*dx1*dx2) * disc

    edx1 = exp(dx1)
    edx2 = exp(dx2)
    
    # 満期タイムステップNにおける原資産価格の初期設定
    S1t = [0 for i in range(2*N+1)] # 2資産あるので2倍の時間(cellに分かれる) 図29の行列
    S2t = [0 for i in range(2*N+1)]

    S1t[0] = S1*exp(-N*dx1) #図29の満期一番下の左下(列でもっとも低下した値)
    S2t[0] = S1*exp(-N*dx2) #図29の満期一番下の左下(行でもっとも低下した値)
    for j in range(1,2*N+1): # j: 1->2N
        S1t[j] = S1t[j-1] * edx1 
        S2t[j] = S2t[j-1] * edx2 

    # 満期における価格のオプション初期設定
    C = np.zeros([N+1, 2*N+1, 2*N+1]) # (2N+1*2N+1)の空行列(3次元配列) 1軸目は時間
    for j in range(0,2*N+1,2): #2ずつ飛ばす 途中はツリーが結合しないので
        for k in range(0,2*N+1,2):
            #C[j][k] = max(0, S1t[j] - S2t[k]) 
            C[N][j][k] = max(0, S1t[j] - S2t[k] - K) 

    # 早期行使条件を適用しながらのツリーのステップバック
    for i in range(0,N)[::-1]: # i:N-1 -> 0 (2,1,0)
        # i=3なら 0,2,4,6, i=2なら1,3,5  i=1なら2,4 i=0なら3
        # -i      -3,-1,1,3       -2,0,2        -1,1       -0 とNずれている
        j = -i+N
        while j <= i+N:
            k = -i+N
            while k <= i+N:
                C[i][j][k] = pdd*C[i+1][j-1][k-1] + pud*C[i+1][j+1][k-1] + pdu*C[i+1][j-1][k+1] + puu*C[i+1][j+1][k+1]
                C[i][j][k] = max(C[i][j][k], S1t[j] - S2t[k] - K)
                k += 2
            j += 2
            
    return C[0][N][N]
