# -*- coding: utf-8 -*-
"""
Created on Tue May 30 16:04:01 2023

@author: nimisha
"""

# LASSOCV
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV
from sklearn.model_selection import TimeSeriesSplit, cross_validate
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
# change the dataset name accordingly
df_final = pd.read_csv('TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall.csv')

df_dataset = df_final.iloc[:, :-1]
df_label = df_final.iloc[:, -1]
cv = RepeatedKFold(n_splits=10, n_repeats=10, random_state=1)


# define model
sel_ = LassoCV(cv=cv, n_jobs=-1)
sel_.fit(df_dataset,df_label)
print('alpha: %f' % sel_.alpha_)



model1 = SelectFromModel(Lasso(alpha=sel_.alpha_,random_state=10))
model = Lasso(alpha=sel_.alpha_,random_state=10)
pipeline = Pipeline(steps=[('s',model1),('m',model)])
model1.fit(df_dataset,df_label)
result = cross_validate(pipeline, df_dataset, df_label, scoring='neg_mean_absolute_error', cv=cv, return_estimator=True)
print(model1.get_support())
print(model1.get_feature_names_out())