import os
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from itertools import combinations
import pickle
from sklearn.preprocessing import PolynomialFeatures
import itertools
import random
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr

# hyper parameter
scale_coef = 100
tol = 1e-6
L0 = 0
L1 = 0.3
L2 = 0.6
H0 = 0.4
H1 = 0.75
H2 = 1.0
ave_approach = "type_ave" # all_ave
# 函数定义
max_min_scalar = lambda x: (x-np.min(x))/(np.max(x)-np.min(x))

# 定义目标函数
def objF1(args):
    Expr_Or,Expr_Rs,Expr_Os = args
    len1 = Expr_Or.shape[1]
    len2 = Expr_Rs.shape[1]
    len3 = Expr_Os.shape[1]
    fx = lambda x: np.mean((x[0]*np.dot(Expr_Or,x[2:2+len1])-x[1]*np.dot(Expr_Rs,x[2+len1:2+len1+len2])-np.dot(Expr_Os,x[2+len1+len2:2+len1+len2+len3]))**2)/scale_coef
    return fx
# constraints1-3 定义函数范围
def make_a_dict_lower(index,Expr_Os,len_Or,len_Rs,lower_lim):
    return {"type":"ineq","fun":lambda x:np.dot(Expr_Os.iloc[index,:],x[2+len_Or+len_Rs:2+len_Or+len_Rs+Expr_Os.shape[1]])- lower_lim[index] }
def make_a_dict_upper(index,Expr_Os,len_Or,len_Rs,upper_lim):
    return {"type":"ineq","fun":lambda x: upper_lim[index] - np.dot(Expr_Os.iloc[index,:],x[2+len_Or+len_Rs:2+len_Or+len_Rs+Expr_Os.shape[1]])}

# constraints 4-10 直接用列表定义

# constraints 11-13 表达量（或平均表达量）的导数为正
def Make_a_derivative_function_Or_fromlist(Ngene,expr,Num_Or,len_Or):
    def make_a_derivative_function_Or_fromlist(x):
        a = x[2:2+Num_Or]
        b = x[2+Num_Or:2+len_Or]
        b_index = 0
        chosen_b = []
        for k in combinations(list(range(Num_Or)),2):
            if Ngene in k:
                chosen_b.append(b[b_index])
            b_index += 1
        chosen_expr = expr[:Num_Or].values.tolist()
        del chosen_expr[Ngene]
        return a[Ngene]+np.dot(chosen_b,chosen_expr)
    return make_a_derivative_function_Or_fromlist
def Make_a_derivative_function_Rs_fromlist(Ngene,expr,Num_Rs,len_Or,len_Rs):
    def make_a_derivative_function_Rs_fromlist(x):
        a = x[2+len_Or:2+len_Or+Num_Rs]
        b = x[2+len_Or+Num_Rs:2+len_Or+len_Rs]
        b_index = 0
        chosen_b = []
        for k in combinations(list(range(Num_Rs)),2):
            if Ngene in k:
                chosen_b.append(b[b_index])
            b_index += 1
        chosen_expr = expr[:Num_Rs].values.tolist()
        del chosen_expr[Ngene]
        return a[Ngene]+np.dot(chosen_b,chosen_expr)
    return make_a_derivative_function_Rs_fromlist
def Make_a_derivative_function_Os_fromlist(Ngene,expr,Num_Os,len_Or,len_Rs):
    def make_a_derivative_function_Os_fromlist(x):
        a = x[2+len_Or+len_Rs:2+len_Or+len_Rs+Num_Os]
        b = x[2+len_Or+len_Rs+Num_Os:]
        b_index = 0
        chosen_b = []
        for k in combinations(list(range(Num_Os)),2):
            if Ngene in k:
                chosen_b.append(b[b_index])
            b_index += 1
        chosen_expr = expr[:Num_Os].values.tolist()
        del chosen_expr[Ngene]
        return a[Ngene]+np.dot(chosen_b,chosen_expr)
    return make_a_derivative_function_Os_fromlist

def truncating(x):
    x[x>np.percentile(x,95)]= np.percentile(x,95)
    return x


# 读入数据
geneExpr = pd.read_csv("all.tpm.16968.1075genes.csv",index_col = 0)
geneExpr = geneExpr*1.0 # DataFrame
geneExpr = geneExpr.apply(truncating,axis = 0) #截尾处理
moduleGene = pd.read_csv("module_genes_input.csv",index_col = 0)
data_gene_all = set(geneExpr.columns)

OrGene = moduleGene.iloc[0,:]
OrGene_ok = []
for i in range(len(OrGene)):
    if OrGene[i] in data_gene_all:
        tmpdata = geneExpr.loc[:,OrGene[i]]
        if (tmpdata<1).sum()/geneExpr.shape[0]<0.3:
            OrGene_ok.append(OrGene[i])
Num_Or = len(OrGene_ok)
RsGene = moduleGene.iloc[1,:]
RsGene_ok = []
for i in range(len(RsGene)):
    if RsGene[i] in data_gene_all:
        tmpdata = geneExpr.loc[:,RsGene[i]]
        if (tmpdata<1).sum()/geneExpr.shape[0]<0.3:
            RsGene_ok.append(RsGene[i])
Num_Rs = len(RsGene_ok)
OsGene = moduleGene.iloc[2,:]
OsGene_ok = []
for i in range(len(OsGene)):
    if OsGene[i] in data_gene_all:
        tmpdata = geneExpr.loc[:,OsGene[i]]
        if (tmpdata<1).sum()/geneExpr.shape[0]<0.3:
            OsGene_ok.append(OsGene[i])
Num_Os = len(OsGene_ok)

Expr_Or = geneExpr[OrGene_ok]
Expr_Rs = geneExpr[RsGene_ok]
Expr_Os = geneExpr[OsGene_ok]
# 带着交互项
poly = PolynomialFeatures(interaction_only=True,include_bias = False) #定义了一个转化器

temp = pd.DataFrame(poly.fit_transform(Expr_Or))
temp.index = Expr_Or.index
original_names = list(Expr_Or.columns)
inter_names = ['_'.join(x) for x in list(itertools.combinations(original_names,2))]
all_names = original_names+inter_names
temp.columns = all_names
Expr_Or0 = temp

temp = pd.DataFrame(poly.fit_transform(Expr_Rs))
temp.index = Expr_Rs.index
original_names = list(Expr_Rs.columns)
inter_names = ['_'.join(x) for x in list(itertools.combinations(original_names,2))]
all_names = original_names+inter_names
temp.columns = all_names
Expr_Rs0 = temp

temp = pd.DataFrame(poly.fit_transform(Expr_Os))
temp.index = Expr_Os.index
original_names = list(Expr_Os.columns)
inter_names = ['_'.join(x) for x in list(itertools.combinations(original_names,2))]
all_names = original_names+inter_names
temp.columns = all_names
Expr_Os0 = temp

Expr_Or0 = Expr_Or0.apply(max_min_scalar,axis = 0)
Expr_Rs0 = Expr_Rs0.apply(max_min_scalar,axis = 0)
Expr_Os0 = Expr_Os0.apply(max_min_scalar,axis = 0)
# 到目前位置 Expr_Or0，Expr_Rs0，Expr_Os0为三个输入矩阵，Num_Or，Num_Rs，Num_Os为三个基因个数，行数为BATCH_SIZE
all_colinfo = pd.read_csv("all.colinfo.csv",index_col = 0)
len_Or = Expr_Or0.shape[1]
len_Rs = Expr_Rs0.shape[1]
len_Os = Expr_Os0.shape[1]


# 以上加载了所有的数据和样本
# 选择样本类型
type_index = [v in ['normal','Solid Tissue Normal','Primary Tumor','Metastatic']for v in all_colinfo.loc[:,"Sample_type"]]
Expr_Or1 = Expr_Or0.loc[type_index,:]
Expr_Rs1 = Expr_Rs0.loc[type_index,:]
Expr_Os1 = Expr_Os0.loc[type_index,:]
colinfo1 = all_colinfo.loc[type_index,:]
labels1  = colinfo1.loc[:,"Sample_type"].values

# 进入循环
coef_array = pd.DataFrame()
all_sample_type = pd.DataFrame()
for epoch in range(100):
    selected_index, _ = train_test_split(list(range(Expr_Or1.shape[0])), test_size=0.88, stratify=labels1)
    Expr_Or = Expr_Or1.iloc[selected_index, :]
    Expr_Rs = Expr_Rs1.iloc[selected_index, :]
    Expr_Os = Expr_Os1.iloc[selected_index, :]
    colinfo = colinfo1.iloc[selected_index, :]
    labels = colinfo.loc[:, "Sample_type"].values
    # 定义F1的上下限
    lower_lim = []
    for i in labels:
        if i == "normal":
            lower_lim.append(L0)
        elif i == "Solid Tissue Normal":
            lower_lim.append(L1)
        else:
            lower_lim.append(L2)
    upper_lim = []
    for i in labels:
        if i == "normal":
            upper_lim.append(H0)
        elif i == "Solid Tissue Normal":
            upper_lim.append(H1)
        else:
            upper_lim.append(H2)
    if ave_approach == "type_ave":
        ave_Or = Expr_Or.groupby(labels).median()
        ave_Rs = Expr_Rs.groupby(labels).median()
        ave_Os = Expr_Os.groupby(labels).median()
    else:
        ave_Or = Expr_Or.mean(axis=0).to_frame().T
        ave_Rs = Expr_Rs.mean(axis=0).to_frame().T
        ave_Os = Expr_Os.mean(axis=0).to_frame().T
    cons_Or_der_ave = []
    for i in range(Num_Or):
        for j in range(ave_Or.shape[0]):
            ave_Or_sub = ave_Or.iloc[j, :]
            cons_Or_der_ave.append(
                {"type": "ineq", "fun": Make_a_derivative_function_Or_fromlist(i, ave_Or_sub, Num_Or, len_Or)})
    cons_Rs_der_ave = []
    for i in range(Num_Rs):
        for j in range(ave_Rs.shape[0]):
            ave_Rs_sub = ave_Rs.iloc[j, :]
            cons_Rs_der_ave.append(
                {"type": "ineq", "fun": Make_a_derivative_function_Rs_fromlist(i, ave_Rs_sub, Num_Rs, len_Or, len_Rs)})
    cons_Os_der_ave = []
    for i in range(Num_Os):
        for j in range(ave_Os.shape[0]):
            ave_Os_sub = ave_Os.iloc[j, :]
            cons_Os_der_ave.append(
                {"type": "ineq", "fun": Make_a_derivative_function_Os_fromlist(i, ave_Os_sub, Num_Os, len_Or, len_Rs)})
    cons_lim = [make_a_dict_lower(i, ave_Os, len_Or, len_Rs, lower_lim) for i in range(ave_Os.shape[0])] + [
        make_a_dict_upper(i, ave_Os, len_Or, len_Rs, upper_lim) for i in range(ave_Os.shape[0])]
    cons_hyper = [{"type": "ineq", "fun": lambda x: 1.5 - x[0]},
                  {"type": "ineq", "fun": lambda x: 1.5 - x[1]},
                  {"type": "ineq", "fun": lambda x: x[0] - 0.5},
                  {"type": "ineq", "fun": lambda x: x[1] - 0.5}]
    x0 = [0.001 for i in range(2 + Expr_Or.shape[1] + Expr_Rs.shape[1] + Expr_Os.shape[1])]
    x0[0] = 1.0
    x0[1] = 1.0
    start = time.time()
    print("begin ........................................")
    res = minimize(objF1((Expr_Or, Expr_Rs, Expr_Os)), x0, method='SLSQP',
                   constraints=cons_lim + cons_hyper + cons_Or_der_ave + cons_Rs_der_ave + cons_Os_der_ave, tol=tol,
                   options={'maxiter': 100, 'ftol': tol, 'disp': True, })
    print(time.time() - start)
    # 进入评估结果阶段
    if epoch == 0:
        coef_array = pd.DataFrame(res.x)
    else:
        coef_array = np.hstack((coef_array, pd.DataFrame(res.x)))
    # 计算所有结果
    values_Os = np.dot(Expr_Os0, res.x[
                                 2 + Expr_Or0.shape[1] + Expr_Rs0.shape[1]:2 + Expr_Or0.shape[1] + Expr_Rs0.shape[1] +
                                                                           Expr_Os0.shape[1]])
    values_Rs = np.dot(Expr_Rs0, res.x[2 + Expr_Or0.shape[1]:2 + Expr_Or0.shape[1] + Expr_Rs0.shape[1]])
    values_Or = np.dot(Expr_Or0, res.x[2:2 + Expr_Or0.shape[1]])
    aaa = np.hstack((values_Or.reshape(Expr_Or0.shape[0], 1), values_Rs.reshape(Expr_Or0.shape[0], 1),
                     values_Os.reshape(Expr_Or0.shape[0], 1),
                     res.x[0] * values_Or.reshape(Expr_Or0.shape[0], 1) - res.x[1] * values_Rs.reshape(
                         Expr_Or0.shape[0], 1)))
    aaa = pd.DataFrame(aaa)
    aaa.columns = ["OR", "RS", "OS", "OOS"]
    aaa.index = all_colinfo.index
    aaa = pd.concat([aaa, all_colinfo], axis=1)
    if epoch == 0:
        all_sample_type = pd.concat((aaa.groupby("Sample_type")["OS"].mean(), aaa.groupby("Sample_type")["RS"].mean(),
                                     aaa.groupby("Sample_type")["OR"].mean()), axis=1)
    else:
        all_sample_type = pd.concat((all_sample_type, aaa.groupby("Sample_type")["OS"].mean(),
                                     aaa.groupby("Sample_type")["RS"].mean(), aaa.groupby("Sample_type")["OR"].mean()),
                                    axis=1)


# 出循环整理结果
coef_array = pd.DataFrame(coef_array)
var = ["alfa","beta"]+Expr_Or.columns.tolist()+Expr_Rs.columns.tolist()+Expr_Os.columns.tolist()
coef_array.index =  var
coef_array.to_csv("./output/coef.csv")
all_sample_type = pd.DataFrame(all_sample_type)
all_sample_type.to_csv("./output/all_sample_type.csv")


