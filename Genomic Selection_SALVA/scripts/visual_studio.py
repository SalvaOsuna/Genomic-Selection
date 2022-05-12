from ctypes import sizeof
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#create an empty 20x20 matrix with row and column labels from 1 to 20:
df = pd.DataFrame(index=range(1, 20), columns=range(1, 20))
print(df)
#In df, when i + j = 20 or i = j print 0s, rest print consecutive numbers from 1 to 324:
for i in range(1, 20):
    for j in range(1, 20):
        if i + j == 20 or i == j:
            df.loc[i, j] = 0
        else:
            df.loc[i, j] = np.random.randint(1, 324,size=324)
print(df)
#create a random 18x18 matrix from 1 to 324:
df_2 = pd.DataFrame(index=range(1, 18), columns=range(1, 18))
for i in range(1, 18):
    for j in range(1, 18):
        df_2.loc[i, j] = np.random.randint(1, 324,size=324)
print(df_2)
