import numpy as np


def get_penalty(x: str, y: str, miss_match: int, plentygap: int):
    i,j,m,n = 0,0,len(x),len(y)
    dp = np.zeros([m + 1, n + 1], dtype=int)
    dp[0:(m + 1), 0] = [i * plentygap for i in range(m + 1)]
    dp[0, 0:(n + 1)] = [i * plentygap for i in range(n + 1)]
    i = 1
    while i <= m:
        j = 1
        while j <= n:
            if x[i - 1] == y[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = min(dp[i - 1][j - 1] + miss_match,dp[i - 1][j] + plentygap,dp[i][j - 1] + plentygap)
            j += 1
        i += 1
    l = n + m
    i = m
    j = n
    posx = l
    posy = l
    xres = np.zeros(l + 1, dtype=int)
    yres = np.zeros(l + 1, dtype=int)
    while not (i == 0 or j == 0):
        if x[i - 1] == y[j - 1]:
            xres[posx] = ord(x[i - 1])
            yres[posy] = ord(y[j - 1])
            posx -= 1
            posy -= 1
            i -= 1
            j -= 1
        elif (dp[i - 1][j - 1] + miss_match) == dp[i][j]:
            xres[posx] = ord(x[i - 1])
            yres[posy] = ord(y[j - 1])
            posx -= 1
            posy -= 1
            i -= 1
            j -= 1
        elif (dp[i - 1][j] + plentygap) == dp[i][j]:
            xres[posx] = ord(x[i - 1])
            yres[posy] = ord('_')
            posx -= 1
            posy -= 1
            i -= 1
        elif (dp[i][j - 1] + plentygap) == dp[i][j]:
            xres[posx] = ord('_')
            yres[posy] = ord(y[j - 1])
            posx -= 1
            posy -= 1
            j -= 1
    while posx > 0:
        if i > 0:
            i -= 1
            xres[posx] = ord(x[i])
            posx -= 1
        else:
            xres[posx] = ord('_')
            posx -= 1
    while posy > 0:
        if j > 0:
            j -= 1
            yres[posy] = ord(y[j])
            posy -= 1
        else:
            yres[posy] = ord('_')
            posy -= 1
    id = 1
    i = l
    while i >= 1:
        if (chr(yres[i]) == '_') and chr(xres[i]) == '_':
            id = i + 1
            break
        i -= 1
    i = id
    x_seq = ""
    while i <= l:
        x_seq += chr(xres[i])
        i += 1
    i = id
    y_seq = ""
    while i <= l:
        y_seq += chr(yres[i])
        i += 1
    print(f"last num in matrix: {dp[m][n]}")
    print(f"X seq: {x_seq}")
    print(f"Y seq: {y_seq}")
    print(f"matrix: "
          f"\n {dp}")

def main():
    gene1 = "AGGGCTTTAAGGACGT"
    gene2 = "AGGCATTTAAGGACG"
    mismatch_penalty = 3
    gap_penalty = 2
    get_penalty(gene1, gene2, mismatch_penalty, gap_penalty)

if __name__ == '__main__':
    main()






