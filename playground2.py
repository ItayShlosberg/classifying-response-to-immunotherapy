import pandas as pd
import pandas
from collections import Counter


# path = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\research\articles\Tables and files article5\Table S1 - A Summary of Data Related to All Single-Cells Analysis, Related to Figure 1.xlsx'
#
#
#
# xls = pd.ExcelFile(path)
# df1 = pd.read_excel(xls, 'Gene marker-Fig1B-C')
#
#
# xls = pd.ExcelFile(path)
# df2 = pd.read_excel(xls, 'Cluster annotation-Fig1B-C')
# print(df2.columns)



x= b.replace('0', '')
x = [x[i:i+20] for i in range(0, len(x)-20, 20)]
print(b[b.index('8')-10: b.index('8')+10])
print(Counter(b.split('\t')))


