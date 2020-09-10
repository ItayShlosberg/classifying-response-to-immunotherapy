import pandas as pd


path = r'C:\Users\itay\Desktop\Technion studies\Keren Laboratory\research\articles\Tables and files article5\Table S1 - A Summary of Data Related to All Single-Cells Analysis, Related to Figure 1.xlsx'



xls = pd.ExcelFile(path)
df1 = pd.read_excel(xls, 'Gene marker-Fig1B-C')


xls = pd.ExcelFile(path)
df2 = pd.read_excel(xls, 'Cluster annotation-Fig1B-C')
print(df2.columns)